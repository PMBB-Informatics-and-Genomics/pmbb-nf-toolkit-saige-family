nextflow.enable.dsl = 2

params.min_survival_cases = 100
params.enable_chunking = false
params.chunks_manifest = null
params.full_chromosome_list = params.chromosome_list
params.host = ""
params.max_vars_for_GRM = null
params.pruning_r2_for_GRM = null
params.min_bin_cases = 50
params.min_quant_n = 500

include { SAIGE_PREPROCESSING } from '../processes/saige_preprocessing.nf'

include { SAIGE_STEP1 } from '../processes/saige_step1.nf'

include { SAIGE_VAR_STEP2 } from '../processes/saige_variant_step2.nf'

include {
    merge_and_filter_saige_gwas_output
    gwas_make_biofilter_positions_input
    make_summary_suggestive_gwas
    make_summary_table_with_annot
    gzipFiles
    MERGE_CHUNKS
    } from '../processes/saige_postprocessing.nf'

include {
    make_pheno_covar_summary_plots
    make_saige_gwas_plots
    make_gwas_plots_with_annot
    collect_gwas_plots
    make_saige_gwas_report
    } from '../processes/saige_visualization.nf'

include { BIOFILTER_POSITIONS } from '../workflows/biofilter_wrapper.nf'

include {
    get_script_file_names
    paramToList
    readChunksManifest
    dump_params_to_json
} from '../processes/saige_helpers.nf'

workflow {
    cohort_pheno_sumstats = SAIGE_GWAS()
}
workflow SAIGE_GWAS {
  main:
    log.info([
        "  NEXTFLOW - DSL2 - SAIGE GWAS - P I P E L I N E",
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
        String.format("  %-25s : %s", "step2_pgen_prefix", params.step2_pgen_prefix),
        "",
        "  Chunking Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "enable_chunking", params.enable_chunking),
        String.format("  %-25s : %s", "full_chromosome_list", params.full_chromosome_list),
        String.format("  %-25s : %s", "chunks_manifest", params.chunks_manifest),
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
        "  SAIGE GWAS Parameters",
        "  " + "=" * 50,
        String.format("  %-25s : %s", "is_Firth_beta", params.use_firth),
        String.format("  %-25s : %s", "pCutoffforFirth", params.firth_cutoff),
        String.format("  %-25s : %s", "LOCO", params.LOCO),
        String.format("  %-25s : %s", "gwas_col_names", params.gwas_col_names),
        String.format("  %-25s : %s", "annotate", params.annotate),
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
        String.format("  %-25s : %s", "my_bgenix", params.my_bgenix),
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

    // subsample chunks down to only those chromosomes specified
    if (params.enable_chunking) {
        assert params.chunks_manifest != null : \
            "enable_chunking is true but chunks_manifest is not set. " +
            "Provide a tab-separated manifest file with columns: chromosome<TAB>file_path_stem"
        def all_chunks = readChunksManifest(params.chunks_manifest.toString())
        filtered_chunks_list = all_chunks.findAll { chrom, stem ->
            params.chromosome_list.any { it.toString() == chrom }
        }
        assert filtered_chunks_list.size() > 0 : \
            "enable_chunking is true but no chunks matched chromosome_list. " +
            "Check that chunks_manifest chromosomes match params.chromosome_list values."
        covered_chroms = filtered_chunks_list.collect { it[0] }.toSet()
        missing_chroms = params.chromosome_list.collect { it.toString() }.findAll { !covered_chroms.contains(it) }
        assert missing_chroms.isEmpty() : \
            "chunks_manifest is missing entries for chromosome(s): ${missing_chroms.join(', ')}. " +
            "Every chromosome in chromosome_list must have at least one chunk in the manifest."
        chunks_list = filtered_chunks_list  // list of [chrom, stem] pairs
    }

    // figure out which samplefile to use
    if (ftype == 'PLINK1') {
        if (params.enable_chunking) {
            step2_fam = "${chunks_list[0][1]}.fam"
        } else {
            step2_fam = "${params.step2_plink_prefix}${params.chromosome_list[0]}.fam"
        }
    } else if (ftype == 'PLINK2') {
        if (params.enable_chunking) {
            step2_fam = "${chunks_list[0][1]}.psam"
        } else {
            step2_fam = "${params.step2_pgen_prefix}${params.chromosome_list[0]}.psam"
        }
    } else if (ftype == 'BGEN') {
        step2_fam = params.bgen_samplefile
    } else {
        //********************Improper File Type****************************
        throw new Exception('Improper file type for step 2, please refer to your .config file \
                            and ensure the ftype is PLINK1, PLINK2, or BGEN. (case matters)')
    }

    // Call Preprocessing sub-workflow (SAIGE_PREPROCESSING)
    workflow_is_phewas = false
    preprocessing_output = SAIGE_PREPROCESSING(pheno_covar_table, cohort_table, step1_fam, step2_fam, workflow_is_phewas)
    keep_cohort_bin_pheno_combos = preprocessing_output[0]
    keep_cohort_quant_pheno_combos = preprocessing_output[1]
    keep_cohort_survival_pheno_combos = preprocessing_output[2]
    pheno_table = preprocessing_output[3]
    cohort_sample_lists = preprocessing_output[4] //sample_list.txt path
    cohort_pheno_tables = preprocessing_output[5] //saige_pheno_covars.txt path

    // Call Step 1 sub-workflow (SAIGE_STEP1)
    step1_is_gene = false
    use_plink_prefix = params.step1_plink_prefix

    (step1_bin_output,step1_quant_output,step1_survival_output) = SAIGE_STEP1(cohort_sample_lists,
          cohort_pheno_tables,
          keep_cohort_bin_pheno_combos,
          keep_cohort_quant_pheno_combos,
          keep_cohort_survival_pheno_combos,
          use_plink_prefix,
          step1_is_gene)

    use_genetic_data_prefix = ftype == 'PLINK1' ? params.step2_plink_prefix \
                            : ftype == 'PLINK2' ? params.step2_pgen_prefix \
                            : params.step2_bgen_prefix
    bgen_sample_file = ftype == 'BGEN' ? params.bgen_samplefile : null
    
    (step2_bin_output, step2_quant_output,step2_survival_output) = SAIGE_VAR_STEP2(
          step1_bin_output,
          step1_quant_output,
          step1_survival_output,
          use_genetic_data_prefix,
          bgen_sample_file,
          workflow_is_phewas
      )
    
    //Add in merging step for chunks if the flag == TRUE
    if (params.enable_chunking) {
        println "Merging chunks"
        step2_all_output = step2_bin_output.concat(step2_quant_output,step2_survival_output)

        //extract the chromosome name from  in tuple
        step2_all_output = step2_all_output.map { cohort_dir, pheno, chr, file_path ->
        def full_chromosome = chr.tokenize('_')[0]
        return [cohort_dir, pheno, full_chromosome, file_path]
        }
      
        //group them by cohort, phenotype, and full chromosome
        step2_grouped_output = step2_all_output.groupTuple(by: [0,1,2])
        // from saige_postprocessing.nf MERGE_CHUNKS process
        step2_grouped_output = MERGE_CHUNKS(step2_grouped_output)
        step2_grouped_output = step2_grouped_output.groupTuple(by: [0, 1], size: params.chromosome_list.size())
        // step2_grouped_output.view{"sgo: ${it}"}

    } else {
        println "No chunks to merge"
        step2_all_output = step2_bin_output.concat(step2_quant_output, step2_survival_output)
        step2_grouped_output = step2_all_output.groupTuple(by: [0, 1], size: params.chromosome_list.size())
    }

    /*
      Step 2 -> Merged Sumstats Channel Emission Tuples
      Step 2 Out:  cohort, phenotype, chromosome, regions, singles
      Group By:    cohort, phenotype
      Merge In:    cohort, phenotype, [chr_list], [region_list], [singles_list]
      - then map to split singles vs regions
    */
    // Collect saige output into channels for merge
    
    merge_singles_script = script_name_dict['merge']

    (singles_merge_output, filtered_singles_output) =  \
    merge_and_filter_saige_gwas_output(step2_grouped_output, merge_singles_script)
    
    filtered_singles_output_list = filtered_singles_output.collect() //make list of filter paths
    pos_input = filtered_singles_output.map { cohort, phecode, filtered_sumstats_path -> new Tuple(filtered_sumstats_path) }.collect()
    
    // ANNOTATIONS
    if (params['annotate']) {
        biofilter_input = gwas_make_biofilter_positions_input(pos_input)
        bf_input_channel = Channel.of('GWAS').combine(biofilter_input)
        biofilter_annots = BIOFILTER_POSITIONS(bf_input_channel)
        plotting_script = script_name_dict['gwas_plots_with_annot']
        anno_input = filtered_singles_output.combine(biofilter_annots)
        plots = make_gwas_plots_with_annot(singles_merge_output.combine(biofilter_annots), \
                                            plotting_script, pheno_table)
        gwas_summary = make_summary_table_with_annot(pos_input, biofilter_annots)
        gwas_csvs = plots.flatten().filter{ it.name.endsWith('manifest.csv') }.collect()
        gwas_analysis = "gwas"
        gwas_manifest = collect_gwas_plots(gwas_analysis, gwas_csvs)
        gwas_pngs = plots.flatten().filter{ it.name =~ /.*png/ }.collect()
    } else {
        gwas_plots_script = script_name_dict['gwas_plots']
        gwas_plots = make_saige_gwas_plots(singles_merge_output, gwas_plots_script, pheno_table)
        gwas_summary = make_summary_suggestive_gwas(pos_input)
        gwas_csvs_2 = gwas_plots.flatten().filter{ it.name.endsWith('manifest.csv') }.collect()
        gwas_analysis_2 = "gwas"
        gwas_manifest = collect_gwas_plots(gwas_analysis_2, gwas_csvs_2)
        gwas_pngs = gwas_plots.flatten().filter{ it.name =~ /.*png/ }.collect()
    }

    json_params = dump_params_to_json(params, 'saige_gwas')

    pheno_pngs = preprocessing_output[6].flatten().collect()

    all_plot_pngs = gwas_pngs
        .concat(pheno_pngs)
        .collect()

    gwas_manifest_script = "${moduleDir}/../scripts/generate_saige_gwas_manifest.py"
    gwas_report_script   = "${moduleDir}/../scripts/generate_saige_gwas_reports.py"
    make_saige_gwas_report(all_plot_pngs, gwas_summary, pheno_table, json_params, gwas_manifest_script, gwas_report_script)

    emit:
    singles_merge_output
    pheno_table
}
