nextflow.enable.dsl = 2

params.ftype = null // This is a GWAS param

log.info """\
    NEXTFLOW - DSL2 - SAIGE ExWAS - P I P E L I N E
    ==================================================
    run as                  : ${workflow.commandLine}
    run location            : ${launchDir}
    started at              : ${workflow.start}
    python exe              : ${params.my_python}

    Cohorts, Phenotypes, and Chromosomes
    ==================================================
    cohort_list             : ${params.cohort_list}
    sex-stratified_cohorts  : ${params.sex_strat_cohort_list}
    bin_pheno_list          : ${params.bin_pheno_list}
    quant_pheno_list        : ${params.quant_pheno_list}
    sex_specific_pheno_file : ${params.sex_specific_pheno_file}
    chromosome_list         : ${params.chromosome_list}
    cat_covars              : ${params.cat_covars}
    cont_covars             : ${params.cont_covars}
    data_csv                : ${params.data_csv}
    cohort_sets             : ${params.cohort_sets}

    Input File Prefixes
    ==================================================
    use_sparse_GRM          : ${params.use_sparse_GRM}
    step1_sparse_grm        : ${params.step1_sparse_grm}
    step1_sparse_grm_samples: ${params.step1_sparse_grm_samples}

    exome_plink_prefix      : ${params.exome_plink_prefix}
    group_file_prefix       : ${params.group_file_prefix}
    gene_location_file      : ${params.gene_location_file}

    SAIGE Step 1 Plink QC Parameters
    ==================================================
    min maf (maf)           : ${params.maf}
    max missingness (geno)  : ${params.geno}
    hardy-weinberg (hwe)    : ${params.hwe}

    SAIGE-GENE Parameters
    ==================================================
    minMAF                  : ${params.min_maf}
    minMAC                  : ${params.min_mac}
    maxMAF_in_groupTest     : ${params.grouptest_maf}
    annotation_in_groupTest : ${params.grouptest_annotation}
    is_Firth_beta           : ${params.use_firth}
    pCutoffforFirth         : ${params.firth_cutoff}
    LOCO                    : ${params.LOCO}

    Other Parameters
    ==================================================
    step1_script            : ${params.step1_script}
    step2_script            : ${params.step2_script}
    pheno_file_id_col       : ${params.id_col}
    p_cutoff_summarize      : ${params.p_cutoff_summarize}

    """.stripIndent()

include { SAIGE_PREPROCESSING } from '../processes/saige_preprocessing.nf'

include { SAIGE_STEP1 } from '../processes/saige_step1.nf'

include { SAIGE_GENE_STEP2 } from '../processes/saige_gene_step2.nf'

include {
    merge_and_filter_saige_gene_regions_output
    merge_and_filter_saige_gene_singles_output
    make_summary_regions_output
    make_summary_singles_output
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_exwas_regions_plots
    make_saige_exwas_singles_plots
    make_exwas_report_src
    make_exwas_report
    make_exwas_report_methods_blurb
    collect_exwas_regions_plots
    collect_exwas_singles_plots
    } from '../processes/saige_visualization.nf'

include {
    paramToList
    get_script_file_names
} from '../processes/saige_helpers.nf'

workflow {
    // Get the script name manifest from the helper functions
    script_name_dict = get_script_file_names()

    // For ExWAS, we use the same merged exome files for step 1 and step 2
    params.step1_plink_prefix = params.exome_plink_prefix

    // Define input file paths
    pheno_covar_table = params.data_csv
    cohort_table = params.cohort_sets
    step1_fam = "${params.step1_plink_prefix}.fam"
    exome_fam = "${params.exome_plink_prefix}.fam"

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = false
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, exome_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    pheno_table = preprocessing_output[2]
    cohort_sample_lists = preprocessing_output[3]
    cohort_pheno_tables = preprocessing_output[4]

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = true
    use_plink_prefix = params.use_sparse_GRM ? params.step1_plink_prefix : params.exome_plink_prefix
    (step1_bin_output, step1_quant_output) = SAIGE_STEP1(cohort_sample_lists,
        cohort_pheno_tables,
        keep_cohort_bin_pheno_combos,
        keep_cohort_quant_pheno_combos,
        use_plink_prefix,
        step1_is_gene)

    // Call Step 2 sub-workflow (SAIGE_GENE_STEP2)
    use_plink_prefix = params.exome_plink_prefix // Use exome prefix for step 2
    (step2_bin_output, step2_quant_output) = SAIGE_GENE_STEP2(
        step1_bin_output,
        step1_quant_output,
        use_plink_prefix,
        workflow_is_phewas)

    /*
    Step 2 -> Merged Sumstats Channel Emission Tuples
    Step 2 Out:  cohort, phenotype, chromosome, regions, singles
    Group By:    cohort, phenotype
    Merge In:    cohort, phenotype, [chr_list], [region_list], [singles_list]
      - then map to split singles vs regions
    */
    // Collect saige output into channels for merge
    step2_all_output = step2_bin_output.concat(step2_quant_output)
    step2_grouped_output = step2_all_output.groupTuple(by: [0, 1], size: params.chromosome_list.size())

    // extract singles results files from the tuples
    singles_sumstats_chr_input = step2_grouped_output.map {
        cohort, pheno, chr, region, singles, marker_list -> \
        new Tuple(cohort, pheno, chr, singles)
    }

    // extract regions results files from the tuples
    regions_sumstats_chr_input = step2_grouped_output.map {
        cohort, pheno, chr, region, singles, marker_list -> \
        new Tuple(cohort, pheno, chr, region)
    }

    merge_script = merge_script = script_name_dict['merge']
    (singles_merge_output, filtered_singles_output, ur_output, filtered_ur_output) = merge_and_filter_saige_gene_singles_output(singles_sumstats_chr_input, merge_script)
    (regions_merge_output, filtered_regions_output, cauchy_output, filtered_cauchy_output) = merge_and_filter_saige_gene_regions_output(regions_sumstats_chr_input, merge_script, pheno_table)

    // collect a list of just the filtered output files, don't need a wildcards anymore
    summary_singles_input = filtered_singles_output.map { cohort, pheno, filtered -> filtered }.collect()
    singles_summary = make_summary_singles_output(summary_singles_input)

    // collect a list of just the filtered output files, don't need a wildcards anymore
    summary_regions_input = filtered_regions_output.map { cohort, pheno, filtered -> filtered }.collect()
    regions_summary = make_summary_regions_output(summary_regions_input)

    /*
    Post-processing:
    Plots (regions and singles)
    Methods blurb with configuration info
    Report file collection (src directory)
    HTML generation
    */

    gene_file = "${params.gene_location_file}"

    // regions plots
    regions_plots_script = script_name_dict['exwas_regions_plots']
    regions_plots = make_saige_exwas_regions_plots(regions_merge_output, gene_file, regions_plots_script, pheno_table)
    // take the 2 input tuples of pngs and csvs, extract csvs, filter on manifest
    regions_csvs = regions_plots.map{pngs, csvs -> new Tuple(csvs)}.transpose().filter{ it.name =~ /.*manifest.csv/ }.collect()
    exwas_regions = "exwas_regions"
    regions_manifest = collect_exwas_regions_plots(exwas_regions, regions_csvs)

    // singles plots
    singles_plots_script = script_name_dict['exwas_singles_plots']
    singles_plots = make_saige_exwas_singles_plots(singles_merge_output, gene_file, singles_plots_script, pheno_table)
    // take the 2 input tuples of pngs and csvs, extract csvs, filter on manifest
    singles_csvs = singles_plots.map{pngs, csvs -> new Tuple(csvs)}.transpose().filter{ it.name =~ /.*manifest.csv/ }.collect()
    exwas_singles = "exwas_singles"
    singles_manifest = collect_exwas_singles_plots(exwas_singles, singles_csvs)

    // methods info
    methods_script = script_name_dict['exwas_methods']
    methods_blurb = make_exwas_report_methods_blurb(methods_script, params)
}
