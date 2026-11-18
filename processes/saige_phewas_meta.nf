// PheWAS-META processes: harmonize each run cohort's SVP exposure-mode sumstats, then meta-analyze per analysis

// test id columns as a list: a command-line value (--test_id_cols exposure,outcome) arrives as a String
List testIdCols() {
    def cols = params.test_id_cols
    return cols instanceof List ? cols.collect { it.toString() } : cols.toString().split(',').collect { it.trim() }
}

// Memory: dynamic, as in saige_postprocessing.nf. Only OOM-type exits (137-140) are retried, with
// base * attempt; any other failure (e.g. a harmonization error, exit 1) stops the run.
// Measured: meta-analysis of 300 exposures x 2000 outcomes x 10 cohorts (6M rows) peaks at ~2 GB.

process harmonize_phewas_sumstats {
    publishDir "${launchDir}/Harmonized/"
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    memory {
        def base_mem = params.host == 'AOU' ? 16.GB : 8.GB
        return base_mem * task.attempt
    }
    input:
        // config_json: the biobank's merged YAML config for this run cohort (see workflows/saige_phewas_meta.nf)
        // sumstats: may be an empty list when every expected file is missing (they are listed in config_json)
        tuple val(meta_cohort), val(config_json), path(sumstats, stageAs: 'sumstats/*'), path(encoding)
        path harmonize_script
    output:
        tuple val(meta_cohort), path("${meta_cohort}.harmonized.tsv.gz"), emit: harmonized
        path "${meta_cohort}.harmonization_log.csv", emit: log
    script:
        """
        mkdir -p sumstats
        cat > config.json <<'JSON_EOF'
${config_json}
JSON_EOF
        ${params.my_python} ${harmonize_script} \\
            --config config.json \\
            --encoding ${encoding} \\
            --sumstats_dir sumstats \\
            --out_prefix ${meta_cohort}
        """
    stub:
        """
        touch ${meta_cohort}.harmonized.tsv.gz ${meta_cohort}.harmonization_log.csv
        """
}

process run_phewas_meta {
    publishDir "${launchDir}/Meta_Sumstats/${analysis}/"
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    memory {
        def base_mem = params.host == 'AOU' ? 16.GB : 8.GB
        return base_mem * task.attempt
    }
    input:
        // cohort_order: the analysis's meta cohort labels, in params.analyses order; harmonized: their long tables
        tuple val(analysis), val(cohort_order), path(harmonized)
        path meta_script
    output:
        tuple val(analysis), path("${analysis}.phewas_meta.tsv.gz"), emit: meta
        tuple val(analysis), path("${analysis}.phewas_meta.effects.tsv.gz"), emit: effects
        tuple val(analysis), path("${analysis}.phewas_meta.plot_input.tsv.gz"), emit: plot_input
        path "${analysis}.phewas_meta.filtered.csv", emit: filtered
        tuple val(analysis), path("${analysis}.phewas_meta.filtered.csv"), emit: filtered_by_analysis
        path "${analysis}.phewas_meta.test_counts.csv", emit: test_counts
        path "${analysis}.phewas_meta.sample_sizes.csv", emit: sample_sizes
    script:
        """
        ${params.my_python} ${meta_script} \\
            --analysis ${analysis} \\
            --cohort_order ${cohort_order.join(' ')} \\
            --inputs ${harmonized} \\
            --pthresh ${params.p_cutoff_summarize} \\
            --min_studies ${params.min_meta_studies} \\
            --p_col ${params.meta_p_col} \\
            --key_cols ${testIdCols().join(' ')}
        """
    stub:
        """
        touch ${analysis}.phewas_meta.tsv.gz ${analysis}.phewas_meta.effects.tsv.gz ${analysis}.phewas_meta.plot_input.tsv.gz
        touch ${analysis}.phewas_meta.filtered.csv ${analysis}.phewas_meta.test_counts.csv ${analysis}.phewas_meta.sample_sizes.csv
        """
}

process plot_phewas_meta {
    publishDir "${launchDir}/Plots/"
    // safe_to_skip: conf/nextflow.config's withLabel errorStrategy (ignore/terminate) replaces this one when loaded
    label 'safe_to_skip'
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : (params.skipPostprocessingErrors ? 'ignore' : 'terminate') }
    memory {
        def base_mem = params.host == 'AOU' ? 16.GB : 8.GB
        return base_mem * task.attempt
    }
    input:
        tuple val(analysis), path(plot_input)
        path plotting_script
        path pheno_descriptions_file
    output:
        path "${analysis}.Plotting.log.txt"
        // optional: an analysis with no tests has no plots
        path "${analysis}.phewas_plot.*.png", optional: true
    script:
        """
        touch Plotting.log.txt
        ${params.my_python} ${plotting_script} \\
            --pheno_descriptions ${pheno_descriptions_file} \\
            --max_labels ${params.varphewas_plot_max_labels} \\
            --max_labels_per_category ${params.varphewas_plot_max_labels_per_category} \\
            --input_file_list ${plot_input}
        mv Plotting.log.txt ${analysis}.Plotting.log.txt
        shopt -s nullglob
        for f in manhattan_plot_variant_*.png; do
            mv "\$f" "${analysis}.phewas_plot.\${f#manhattan_plot_variant_}"
        done
        """
    stub:
        """
        touch ${analysis}.Plotting.log.txt ${analysis}.phewas_plot.stub.png
        """
}

process forest_plots_phewas_meta {
    publishDir "${launchDir}/Plots/"
    // safe_to_skip: conf/nextflow.config's withLabel errorStrategy (ignore/terminate) replaces this one when loaded
    label 'safe_to_skip'
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : (params.skipPostprocessingErrors ? 'ignore' : 'terminate') }
    memory {
        def base_mem = params.host == 'AOU' ? 8.GB : 4.GB
        return base_mem * task.attempt
    }
    input:
        tuple val(analysis), val(cohort_order), path(filtered), path(effects)
        path forest_script
    output:
        // optional: an empty filtered file (e.g. too few cohorts) gives no plots
        path "${analysis}.forest.*.png", optional: true
    script:
        """
        ${params.my_python} ${forest_script} \\
            --analysis ${analysis} \\
            --filtered ${filtered} \\
            --effects ${effects} \\
            --cohort_order ${cohort_order.join(' ')} \\
            --max_plots ${params.forest_plot_max}
        """
    stub:
        """
        touch ${analysis}.forest.stub.png
        """
}

process make_phewas_meta_report {
    publishDir "${launchDir}", mode: 'copy'
    input:
        path all_plots, stageAs: 'Plots/*'
        path suggestive
        path test_counts
        path sample_sizes
        path harmonization_log
        path cohorts
        path params_json
        path manifest_script
        path report_script
    output:
        path('PheWAS_Meta_Report.zip')
    script:
        """
        ${params.my_python} ${manifest_script} \\
            --params_json ${params_json} \\
            --suggestive ${suggestive} \\
            --test_counts ${test_counts} \\
            --sample_sizes ${sample_sizes} \\
            --harmonization_log ${harmonization_log} \\
            --cohorts ${cohorts} \\
            --plots_dir Plots/ \\
            --output phewas_meta_manifest.json

        ${params.my_python} ${report_script} \\
            --manifest phewas_meta_manifest.json \\
            --output_zip PheWAS_Meta_Report.zip
        """
    stub:
        """
        touch PheWAS_Meta_Report.zip
        """
}
