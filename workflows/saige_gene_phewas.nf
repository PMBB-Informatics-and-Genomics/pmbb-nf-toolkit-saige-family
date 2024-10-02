nextflow.enable.dsl = 2

MIN_BIN_CASES = 50
MIN_QUANT_N = 500

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
    gene_list_file          : ${params.gene_list_file}
    chromosome_list         : ${params.chromosome_list}
    cat_covars              : ${params.cat_covars}
    cont_covars             : ${params.cont_covars}
    data_csv                : ${params.data_csv}
    cohort_sets             : ${params.cohort_sets}
    sex_specific_pheno_file : ${params.sex_specific_pheno_file}

    Input File Prefixes
    ==================================================
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
    merge_and_filter_saige_gene_regions_phewas_output
    merge_and_filter_saige_gene_singles_phewas_output
    make_summary_regions_output
    make_summary_singles_output
    } from '../processes/saige_postprocessing.nf'

include {
    make_saige_gene_phewas_regions_plots
    collect_gene_phewas_regions_plots
} from '../processes/saige_visualization.nf'

include {
    paramToList
    get_script_file_names
} from '../processes/saige_helpers.nf'

workflow {
    // Get the script name manifest from the helper functions
    script_name_dict = get_script_file_names()

    // For Gene PheWAS, we use the same merged exome files for step 1 and step 2
    params.step1_plink_prefix = params.exome_plink_prefix

    pheno_covar_table = params.data_csv
    cohort_table = params.cohort_sets
    step1_fam = "${params.step1_plink_prefix}.fam"
    exome_fam = "${params.exome_plink_prefix}.fam"

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = true
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, exome_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    pheno_table = preprocessing_output[2]
    cohort_sample_lists = preprocessing_output[3]
    cohort_pheno_tables = preprocessing_output[4]

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = true
    use_plink_prefix = params.exome_plink_prefix
    (step1_bin_output, step1_quant_output) = SAIGE_STEP1(cohort_sample_lists,
        cohort_pheno_tables,
        keep_cohort_bin_pheno_combos,
        keep_cohort_quant_pheno_combos,
        use_plink_prefix,
        step1_is_gene)

    // Call Step 2 sub-workflow (SAIGE_GENE_STEP2)
    use_plink_prefix = params.exome_plink_prefix
    (step2_bin_output, step2_quant_output) = SAIGE_GENE_STEP2(
        step1_bin_output,
        step1_quant_output,
        use_plink_prefix,
        workflow_is_phewas)

    /*
    Step 2 -> Merged Sumstats Channel Emission Tuples
    Step 2 Out:  cohort, phenotype, chromosome, regions, singles, marker_list
    Group By:    cohort
    Merge In:    cohort, [phenotype_list], [chr_list], [region_list], [singles_list]
    - then map to split singles vs regions
    */

    // Collect saige output into channels for merge
    step2_all_output = step2_bin_output.concat(step2_quant_output)
    step2_grouped_output = step2_all_output.groupTuple(by: [0, 2])

    // extract singles results files from the tuple
    singles_sumstats_input = step2_grouped_output.map {
        cohort, pheno_list, chr, region, singles, marker_list -> \
        new Tuple(cohort, chr, singles)
    }

    // extract regions results files from the tuple
    regions_sumstats_input = step2_grouped_output.map {
        cohort, pheno_list, chr, region, singles, marker_list -> \
        new Tuple(cohort, chr, region)
    }

    merge_script = script_name_dict['merge']
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_gene_singles_phewas_output(singles_sumstats_input, merge_script)
    (regions_merge_output, filtered_regions_output) = merge_and_filter_saige_gene_regions_phewas_output(regions_sumstats_input, merge_script, pheno_table)

    summary_singles_input = filtered_singles_output.map { cohort, pheno, filtered -> filtered }.collect()
    singles_summary = make_summary_singles_output(summary_singles_input)

    summary_regions_input = filtered_regions_output.map { cohort, pheno, filtered -> filtered }.collect()
    regions_summary = make_summary_regions_output(summary_regions_input)

    // Making Plots
    gene_file = params.gene_location_file
    phenotype_description_file = params.pheno_descriptions_file
    regions_plot_script = "${launchDir}/scripts/make_saige_gene_burden_phewas_plots_v2.py"
    regions_plots = make_saige_gene_phewas_regions_plots(regions_merge_output, phenotype_description_file, regions_plot_script)
    // take the 2 input tuples of pngs and csvs, extract csvs, filter on manifest
    regions_csvs = regions_plots.map{pngs, csvs -> new Tuple(csvs)}.transpose().filter{ it.name =~ /.*manifest.csv/ }.collect()
    phewas_regions = "phewas_regions"
    regions_manifest = collect_gene_phewas_regions_plots(phewas_regions, regions_csvs)
}
