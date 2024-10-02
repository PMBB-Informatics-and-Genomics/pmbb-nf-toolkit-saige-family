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
    GPU                     : ${params.GPU}
    cohort_sets             : ${params.cohort_sets}

    Input File Prefixes
    ==================================================
    step1_plink_prefix      : ${params.step1_plink_prefix}
    step1_sparse_grm        : ${params.step1_sparse_grm}
    step1_sparse_grm_samples: ${params.step1_sparse_grm_samples}
    sparse_grm              : ${params.use_sparse_GRM}
    step2_plink_prefix      : ${params.plink_prefix}
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

include { SAIGE_VAR_STEP2 } from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gene_regions_output
    merge_and_filter_saige_gene_singles_output
    merge_and_filter_saige_gene_singles_phewas_output
    merge_and_filter_saige_variantphewas_output
    make_summary_regions_output
    make_summary_singles_output
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_variantphewas_plots
    } from '../processes/saige_visualization.nf'

include {
    get_script_file_names
} from '../processes/saige_helpers.nf'

workflow {
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
        step2_fam = "${params.step2_plink_prefix}.fam"
    } else if (ftype == 'BGEN') {
        step2_fam = params.samplefile
    } else {
        //********************Improper File Type****************************
        throw new Exception("Improper file type for step 2, please refer to your .config file \
                           and ensure the ftype is PLINK, or BGEN.")
    }

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = true
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, step2_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    pheno_table = preprocessing_output[2]
    cohort_sample_lists = preprocessing_output[3]
    cohort_pheno_tables = preprocessing_output[4]

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = false
    use_plink_prefix = params.step1_plink_prefix
    (step1_bin_output, step1_quant_output) = SAIGE_STEP1(cohort_sample_lists,
        cohort_pheno_tables,
        keep_cohort_bin_pheno_combos,
        keep_cohort_quant_pheno_combos,
        use_plink_prefix,
        step1_is_gene)

    use_genetic_data_prefix = ftype == "PLINK" ? params.step2_plink_prefix : params.step2_bgen_prefix
    bgen_sample_file = ftype == "PLINK" ? null : params.samplefile
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
    Merge In:    cohort,chr_list, [phenotype], [region_list], [singles_list]
      - then map to split singles vs regions
    */
    // Collect saige output into channels for merge
    step2_all_output = step2_bin_output.concat(step2_quant_output)
    step2_grouped_output = step2_all_output.groupTuple(by: [0,2])
    merge_singles_script = "${launchDir}/scripts/merge_and_filter_saige_results.py"
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_variantphewas_output(step2_grouped_output, merge_singles_script)

    // collect a list of just the filtered output files, don't need a wildcards anymore
    summary_singles_input = filtered_singles_output.map { cohort, pheno, filtered -> filtered }.collect()
    singles_summary = make_summary_singles_output(summary_singles_input)
    
    if (params.pheno_descriptions_file) {
        single_varphewas_plot_script = "${launchDir}/scripts/make_saige_var_phewas_plots.py"
        singles_plot_input = singles_merge_output.map { cohort,pheno,path -> path }.collect()
        singles_plots = make_saige_variantphewas_plots(singles_plot_input,single_varphewas_plot_script,params.pheno_descriptions_file) }
    else {
        println("No Phenotype Information given. Plots not generated.")
    }
}
