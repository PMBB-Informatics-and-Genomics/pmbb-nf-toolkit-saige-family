params.host = ""
params.max_vars_for_GRM = null
params.pruning_r2_for_GRM = null
params.min_bin_cases = 50
params.min_quant_n = 500
params.min_survival_cases = 50

include { SAIGE_PREPROCESSING } from '../processes/saige_preprocessing.nf'

include { SAIGE_STEP1 } from '../processes/saige_step1.nf'

include { SAIGE_VAR_STEP2 } from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gene_regions_output
    merge_and_filter_saige_gene_singles_output
    merge_and_filter_saige_gene_singles_phewas_output
    merge_and_filter_saige_variant_phewas_output
    make_summary_regions_output
    make_summary_singles_output
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_variant_phewas_plots
    make_saige_variant_phewas_report
    } from '../processes/saige_visualization.nf'

include {
    get_script_file_names
    dump_params_to_json
} from '../processes/saige_helpers.nf'

workflow {
    log.info([
        "  NEXTFLOW - DSL2 - SAIGE Variant PheWAS - P I P E L I N E",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "run as", workflow.commandLine),
        String.format("  %-25s : %s", "run location", launchDir),
        String.format("  %-25s : %s", "started at", workflow.start),
        String.format("  %-25s : %s", "python exe", params.my_python),
        "",
        "  Cohorts, Phenotypes, and Chromosomes",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "cohort_list", params.cohort_list),
        String.format("  %-25s : %s", "sex-stratified_cohorts", params.sex_strat_cohort_list),
        String.format("  %-25s : %s", "bin_pheno_list", params.bin_pheno_list),
        String.format("  %-25s : %s", "quant_pheno_list", params.quant_pheno_list),
        String.format("  %-25s : %s", "survival_pheno_list", params.survival_pheno_list),
        String.format("  %-25s : %s", "sex_specific_pheno_file", params.sex_specific_pheno_file),
        String.format("  %-25s : %s", "chromosome_list", params.chromosome_list),
        String.format("  %-25s : %s", "cat_covars", params.cat_covars),
        String.format("  %-25s : %s", "cont_covars", params.cont_covars),
        String.format("  %-25s : %s", "sex_strat_cat_covars", params.sex_strat_cat_covars),
        String.format("  %-25s : %s", "sex_strat_cont_covars", params.sex_strat_cont_covars),
        String.format("  %-25s : %s", "data_csv", params.data_csv),
        String.format("  %-25s : %s", "cohort_sets", params.cohort_sets),
        String.format("  %-25s : %s", "min_bin_cases", params.min_bin_cases),
        String.format("  %-25s : %s", "min_quant_n", params.min_quant_n),
        String.format("  %-25s : %s", "min_survival_cases", params.min_survival_cases),
        "",
        "  Input File Prefixes",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "ftype", params.ftype),
        String.format("  %-25s : %s", "use_sparse_GRM", params.use_sparse_GRM),
        String.format("  %-25s : %s", "step1_sparse_grm", params.step1_sparse_grm),
        String.format("  %-25s : %s", "step1_sparse_grm_samples", params.step1_sparse_grm_samples),
        String.format("  %-25s : %s", "step1_plink_prefix", params.step1_plink_prefix),
        String.format("  %-25s : %s", "step2_plink_prefix", params.step2_plink_prefix),
        String.format("  %-25s : %s", "step2_bgen_prefix", params.step2_bgen_prefix),
        String.format("  %-25s : %s", "bgen_samplefile", params.bgen_samplefile),
        String.format("  %-25s : %s", "pheno_descriptions_file", params.pheno_descriptions_file),
        "",
        "  SAIGE Step 1 Plink QC Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "min maf (maf)", params.maf),
        String.format("  %-25s : %s", "max missingness (geno)", params.geno),
        String.format("  %-25s : %s", "hardy-weinberg (hwe)", params.hwe),
        String.format("  %-25s : %s", "max_vars_for_GRM", params.max_vars_for_GRM),
        String.format("  %-25s : %s", "min_vars_for_GRM", params.min_vars_for_GRM),
        String.format("  %-25s : %s", "min_rare_vars_for_GRM", params.min_rare_vars_for_GRM),
        String.format("  %-25s : %s", "pruning_r2_for_GRM", params.pruning_r2_for_GRM),
        "",
        "  SAIGE Variant Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "is_Firth_beta", params.use_firth),
        String.format("  %-25s : %s", "pCutoffforFirth", params.firth_cutoff),
        String.format("  %-25s : %s", "LOCO", params.LOCO),
        String.format("  %-25s : %s", "singles_col_names", params.singles_col_names),
        "",
        "  Survival Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "event_time_col", params.event_time_col),
        String.format("  %-25s : %s", "event_time_bin", params.event_time_bin),
        "",
        "  Other Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "host", params.host),
        String.format("  %-25s : %s", "GPU", params.GPU),
        String.format("  %-25s : %s", "step1_script", params.step1_script),
        String.format("  %-25s : %s", "step2_script", params.step2_script),
        String.format("  %-25s : %s", "pheno_file_id_col", params.id_col),
        String.format("  %-25s : %s", "p_cutoff_summarize", params.p_cutoff_summarize),
    ].join("\n"))

    // cohort = Channel.fromList(params.cohort_list)
    chromosome = Channel.fromList(params.chromosome_list)
    plink_suffixes_list = ['.bed', '.bim', '.fam']
    ftype = params.ftype

    // Get the script name manifest from the helper functions
    script_name_dict = get_script_file_names()

    pheno_covar_table = params.data_csv
    cohort_table = params.cohort_sets
    step1_fam = "${params.step1_plink_prefix}.fam"
    
    if (ftype == 'PLINK1') {
        step2_fam = "${params.step2_plink_prefix}.fam"
    } else if (ftype == 'BGEN') {
        step2_fam = params.bgen_samplefile
    } else {
        //********************Improper File Type****************************
        throw new Exception("Improper file type for step 2, please refer to your .config file \
                           and ensure the ftype is PLINK1, or BGEN.")
    }

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = true
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, step2_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    keep_cohort_survival_pheno_combos = preprocessing_output[2]
    pheno_table = preprocessing_output[3]
    cohort_sample_lists = preprocessing_output[4]
    cohort_pheno_tables = preprocessing_output[5]

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = false
    use_plink_prefix = params.step1_plink_prefix
    (step1_bin_output, step1_quant_output,step1_survival_output) = SAIGE_STEP1(cohort_sample_lists,
        cohort_pheno_tables,
        keep_cohort_bin_pheno_combos,
        keep_cohort_quant_pheno_combos,
        keep_cohort_survival_pheno_combos,
        use_plink_prefix,
        step1_is_gene)

    use_genetic_data_prefix = ftype == 'PLINK1' ? params.step2_plink_prefix : params.step2_bgen_prefix
    bgen_sample_file = ftype == 'BGEN' ? params.bgen_samplefile : null
    (step2_bin_output, step2_quant_output) = SAIGE_VAR_STEP2(
        step1_bin_output,
        step1_quant_output,
        step1_survival_output,
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
    step2_grouped_output = step2_all_output.groupTuple(by: [0, 2])
    merge_singles_script = script_name_dict['merge']
    (singles_merge_output, filtered_singles_output) = merge_and_filter_saige_variant_phewas_output(step2_grouped_output, merge_singles_script)

    // collect a list of just the filtered output files, don't need a wildcards anymore
    summary_singles_input = filtered_singles_output.map { cohort, pheno, filtered -> filtered }.collect()
    singles_summary = make_summary_singles_output(summary_singles_input)

    json_params = dump_params_to_json(params, 'saige_variant_phewas')

    var_phewas_manifest_script = "${moduleDir}/../scripts/generate_saige_variant_phewas_manifest.py"
    var_phewas_report_script   = "${moduleDir}/../scripts/generate_saige_variant_phewas_reports.py"

    if (params.pheno_descriptions_file) {
        single_varphewas_plot_script = "${moduleDir}/../scripts/make_saige_var_phewas_plots.py"
        singles_plot_input = singles_merge_output.groupTuple(by: 0)
        (plot_logs, variant_pngs) = make_saige_variant_phewas_plots(
            singles_plot_input,
            single_varphewas_plot_script,
            params.pheno_descriptions_file
        )
        all_plot_pngs = variant_pngs.collect()
        make_saige_variant_phewas_report(
            all_plot_pngs,
            singles_summary,
            pheno_table,
            json_params,
            var_phewas_manifest_script,
            var_phewas_report_script
        )
    } else {
        println('No Phenotype Information given. Plots not generated. Report not generated.')
    }
}
