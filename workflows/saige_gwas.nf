nextflow.enable.dsl = 2

MIN_BIN_CASES = 50
MIN_QUANT_N = 500

/*
Chris Carson modifying and adding to the work of
Florian Wuennemann and Lindsay Guare for LPC at Penn
*/

log.info """\
    NEXTFLOW - DSL2 - SAIGE GWAS - P I P E L I N E
    ==================================================
    run as                  : ${workflow.commandLine}
    run location            : ${launchDir}
    started at              : ${workflow.start}
    python exe              : ${params.my_python}
    bgenix exe              : ${params.my_bgenix}
    filetype(bgen,plink)    : ${params.ftype}

    Cohorts, Phenotypes, and Chromosomes
    ==================================================
    cohort_list             : ${params.cohort_list}
    bin_pheno_list          : ${params.bin_pheno_list}
    quant_pheno_list        : ${params.quant_pheno_list}
    chromosome_list         : ${params.chromosome_list}
    cat_covars              : ${params.cat_covars}
    cont_covars             : ${params.cont_covars}

    Input File Prefixes
    ==================================================
    use_sparse_GRM          : ${params.use_sparse_GRM}
    step1_plink_prefix      : ${params.step1_plink_prefix}
    step2_plink_prefix     : ${params.step2_plink_prefix}
    step2_bgen_prefix      : ${params.step2_bgen_prefix}

    SAIGE Step 1 Plink QC Parameters
    ==================================================
    min maf (maf)           : ${params.maf}
    max missingness (geno)  : ${params.geno}
    hardy-weinberg (hwe)    : ${params.hwe}

    SAIGE-GENE Parameters
    ==================================================
    minMAF                  : ${params.min_maf}
    minMAC                  : ${params.min_mac}
    pCutoffforFirth         : ${params.firth_cutoff}

    Other Parameters
    ==================================================
    step1_script            : ${params.step1_script}
    step2_script            : ${params.step2_script}
    bgen_samplefile         : ${params.bgen_samplefile}
    pheno_file_id_col       : ${params.id_col}
    p_cutoff_summarize      : ${params.p_cutoff_summarize}
    gwas_col_names          : ${params.gwas_col_names}
    annotate                : ${params.annotate}
    biofilter_loki          : ${params.biofilter_loki}
    biofilter_script        : ${params.biofilter_script}

    """.stripIndent()

include { SAIGE_PREPROCESSING } from '../processes/saige_preprocessing.nf'

include { SAIGE_STEP1 } from '../processes/saige_step1.nf'

include { SAIGE_VAR_STEP2 } from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gwas_output
    gwas_make_biofilter_positions_input
    make_summary_suggestive_gwas
    make_summary_table_with_annot
    gzipFiles
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_gwas_plots
    make_gwas_plots_with_annot
    collect_gwas_plots
    } from '../processes/saige_visualization.nf'

if (params.annotate) {
  include { BIOFILTER_POSITIONS } from "${launchDir}/biofilter_wrapper.nf"
}

include {
  get_script_file_names
} from '../processes/saige_helpers.nf'

workflow {
    cohort_pheno_sumstats = SAIGE_GWAS()
}

workflow SAIGE_GWAS {
  main:
    // cohort = Channel.fromList(params.cohort_list)
    chromosome = Channel.fromList(params.chromosome_list)
    plink_suffixes_list = ['.bed', '.bim', '.fam']
    ftype = params.ftype

    // Get the script name manifest from the helper functions
    script_name_dict = get_script_file_names()

    pheno_covar_table = params.data_csv
    cohort_table = params.cohort_sets
    step1_fam = "${params.step1_plink_prefix}.fam"

    if (ftype == 'PLINK') {
        step2_fam = "${params.step2_plink_prefix}${params.chromosome_list[0]}.fam"
      } else if (ftype == 'BGEN') {
        step2_fam = params.bgen_samplefile
      } else {
        //********************Improper File Type****************************
        throw new Exception('Improper file type for step 2, please refer to your .config file \
                            and ensure the ftype is PLINK, or BGEN. (case matters)')
    }

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = false
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, step2_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    pheno_table = preprocessing_output[2]
    cohort_sample_lists = preprocessing_output[3] //sample_list.txt path
    cohort_pheno_tables = preprocessing_output[4] //saige_pheno_covars.txt path

    // keep_cohort_bin_pheno_combos.view{"kcbp: ${it}"}
    // pheno_table.view{"pt: ${it}"}
    // cohort_pheno_tables.view{"cpt: ${it}"}

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = false
    use_plink_prefix = params.step1_plink_prefix

    (step1_bin_output, step1_quant_output) = SAIGE_STEP1(cohort_sample_lists,
          cohort_pheno_tables,
          keep_cohort_bin_pheno_combos,
          keep_cohort_quant_pheno_combos,
          use_plink_prefix,
          step1_is_gene)

    use_genetic_data_prefix = ftype == 'PLINK' ? params.step2_plink_prefix : params.step2_bgen_prefix
    bgen_sample_file = ftype == 'PLINK' ? null : params.bgen_samplefile
    
    (step2_bin_output, step2_quant_output) = SAIGE_VAR_STEP2(
          step1_bin_output,
          step1_quant_output,
          use_genetic_data_prefix,
          bgen_sample_file,
          workflow_is_phewas
      )
    
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
    merge_singles_script = script_name_dict['merge']

    //step2outputfilestozip = step2_grouped_output.map{cohort_dir, pheno, chr_list, chr_inputs -> new Tuple(chr_inputs.join(' '))}
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_gwas_output(step2_grouped_output, merge_singles_script)
    //sto = chr_inputs.map{chr_inputs-> new Tuple(chr_inputs.join(' '))}
    //gzipFiles(sto)

    filtered_singles_output_list = filtered_singles_output.collect() //make list of filter paths
    pos_input = filtered_singles_output.map { cohort, phecode, filtered_sumstats_path -> new Tuple(filtered_sumstats_path) }.collect()
    // ANNOTATIONS
    if (params['annotate']) {
            biofilter_input = gwas_make_biofilter_positions_input(pos_input)
            bf_input_channel = Channel.of('GWAS').combine(biofilter_input)
            biofilter_annots = BIOFILTER_POSITIONS(bf_input_channel)
            plotting_script = script_name_dict['gwas_plots_with_annot']
            anno_input = filtered_singles_output.combine(biofilter_annots)
            plots = make_gwas_plots_with_annot(singles_merge_output.combine(biofilter_annots), plotting_script, pheno_table)
            make_summary_table_with_annot(pos_input, biofilter_annots)
            gwas_csvs = plots.filter{ it.name =~ /.*manifest.csv/ }.collect()
            gwas_analysis = "gwas"
            // gwas_manifest = collect_gwas_plots(gwas_analysis, gwas_csvs)
    }

    else {
        gwas_plots_script = script_name_dict['gwas_plots']
        gwas_plots = make_saige_gwas_plots(singles_merge_output, gwas_plots_script, pheno_table)
        // filtered_singles_output_list.view { "filtered: ${it }"}
        make_summary_suggestive_gwas(pos_input)
        gwas_csvs_2 = gwas_plots.filter{ it.name =~ /.*manifest.csv/ }.collect()
        gwas_analysis_2 = "gwas"
        // gwas_manifest_2 = collect_gwas_plots(gwas_analysis_2, gwas_csvs_2)
     }

    emit:
      singles_merge_output
      pheno_table
}
