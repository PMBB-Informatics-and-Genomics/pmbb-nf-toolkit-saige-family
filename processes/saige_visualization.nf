import groovy.json.JsonBuilder

process make_pheno_covar_summary_plots {
    publishDir "${launchDir}/Plots/"
    input:
        val cohort_list
        val bin_pheno_list
        val quant_pheno_list
        val survival_pheno_list
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(exome_fam, stageAs: 'Exome/*')
        path pheno_covar_table
        path cohort_table

        path(pheno_covar_plots_script)
    output:
        path('*.png')
    shell:
        """
        ${params.my_python} ${pheno_covar_plots_script} \
          -c ${cohort_list.join(' ')} \
          -b ${bin_pheno_list.join(' ')} \
          -q ${quant_pheno_list.join(' ')} \
          -t ${survival_pheno_list.join(' ')} \
          --step1Fam ${step1_fam} \
          --exomeFam ${exome_fam} \
          --data ${pheno_covar_table} \
          --samples ${cohort_table} \
          --id ${params.id_col}
        """
    stub:
        '''
        touch stub_plot.png
        '''
}

process make_saige_gwas_plots {
    publishDir "${launchDir}/Plots/"
    label 'safe_to_skip', 'high_memory_plots'
    // needs dynamic memory {} allocation
    // Peak memory ~ compressed_size * 10 (uncompressed) + 4GB overhead for plotting.
    memory {
        def fileSizeGb = sumstats.size() / (1024**3)
        def max_mem_gb = 256
        def requiredMemGb = Math.max(16, (Math.ceil(fileSizeGb * 20) + 4).toLong())
        def attempt_mem = Math.min(requiredMemGb, max_mem_gb)
        return attempt_mem.GB
    }

    input:
        // from singles_merge_output
        tuple val(cohort), val(pheno), path(sumstats)
        //
        path(gwas_plot_script)

        path(pheno_table)
    output:
        path "${cohort}.${pheno}.{manhattan_vertical.png,qq.png,qq.csv,gwas.plots_manifest.csv}"
    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${gwas_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -s ${sumstats} \
          -t ${pheno_table}
        """
    stub:
        """
        touch ${cohort}.${pheno}.singles.manhattan_vertical.png
        touch ${cohort}.${pheno}.singles.qq.png
        touch ${cohort}.${pheno}.gwas.plots_manifest.csv
        """
}

process make_gwas_plots_with_annot {
    publishDir "${launchDir}/Plots/"
    memory {
        def fileSizeGb = sumstats.size() / (1024**3)
        def max_mem_gb = 256
        def requiredMemGb = Math.max(16, (Math.ceil(fileSizeGb * 20) + 4).toLong())
        def attempt_mem = Math.min(requiredMemGb, max_mem_gb)
        return attempt_mem.GB
    }
    label 'safe_to_skip', 'high_memory_plots'
    input:
        tuple val(cohort), val(pheno), path(sumstats), val(analysis), path(biofilter_annots)
        path plotting_script
        path pheno_table
    output:
        path "${cohort}.${pheno}.{manhattan_vertical.png,qq.png,qq.csv,gwas.plots_manifest.csv}"
    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${plotting_script} \
          --phenotype ${pheno} \
          --analysis ${analysis} \
          --sumstats ${sumstats} \
          --annot ${biofilter_annots} \
          --cohort ${cohort} \
          --phenoTable ${pheno_table}
        """
    stub:
        """
        touch ${pheno}.${analysis}.manhattan.png
        touch ${pheno}.${analysis}.qq.png
        touch ${pheno}.${analysis}.qq.csv
        touch ${cohort}.${pheno}.gwas.plots_manifest.csv
        """
}

process make_saige_variant_phewas_plots {
    publishDir "${launchDir}/Plots/"
    memory {
        def fileSizeGb = singles_outputs.size() / (1024**3)
        def max_mem_gb = 256
        def requiredMemGb = Math.max(16, (Math.ceil(fileSizeGb * 20) + 4).toLong())
        def attempt_mem = Math.min(requiredMemGb, max_mem_gb)
        return attempt_mem.GB
    }
    label 'safe_to_skip', 'high_memory_plots'
    input:
        tuple val(cohort), val(pheno_list), path(singles_outputs)
        path plotting_script
        path pheno_descriptions_file
    output:
        path 'Plotting.log.txt'
        path '*.manhattan_plot_variant_*.png'
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        echo "${singles_outputs}"
        ${params.my_python} ${plotting_script} \
          --pheno_descriptions ${pheno_descriptions_file} \
          --input_file_list ${singles_outputs}
        for f in manhattan_plot_variant_*.png; do
            mv "\$f" "${cohort}.\$f"
        done
        """
    stub:
        """
        touch Plotting.log.txt
        touch ${cohort}.manhattan_plot_variant_stub.png
        """
}

process make_saige_exwas_singles_plots {
    publishDir "${launchDir}/Plots/"
    memory {
        def fileSizeGb = singles_sumstats.size() / (1024**3)
        def max_mem_gb = 256
        def requiredMemGb = Math.max(16, (Math.ceil(fileSizeGb * 20) + 4).toLong())
        def attempt_mem = Math.min(requiredMemGb, max_mem_gb)
        return attempt_mem.GB
    }
    label 'safe_to_skip', 'high_memory_plots'
    input:
        tuple val(cohort), val(pheno), path(singles_sumstats)

        path(gene_file)
        path(exwas_singles_plot_script)

        path(pheno_table)
    output:
        tuple path("${cohort}.${pheno}*.png"), path("${cohort}.${pheno}*.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${exwas_singles_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -g ${gene_file} \
          -s ${singles_sumstats} \
          -t ${pheno_table} \
          -m ${params.grouptest_maf} \
          -a "${params.grouptest_annotation}"
        """
    stub:
        """
        touch ${cohort}.${pheno}.singles.manhattan_vertical.png
        touch ${cohort}.${pheno}.singles.qq.png
        touch ${cohort}.${pheno}.exwas_singles.plots_manifest.csv
        """
}

process make_saige_exwas_regions_plots {
    publishDir "${launchDir}/Plots/"
    memory {
        def fileSizeGb = regions_sumstats.size() / (1024**3)
        def max_mem_gb = 256
        def requiredMemGb = Math.max(16, (Math.ceil(fileSizeGb * 20) + 4).toLong())
        def attempt_mem = Math.min(requiredMemGb, max_mem_gb)
        return attempt_mem.GB
    }
    label 'safe_to_skip', 'high_memory_plots'
    input:
        tuple val(cohort), val(pheno), path(regions_sumstats)

        path(gene_file)
        path(exwas_regions_plot_script)

        path(pheno_table)
    output:
        tuple path("${cohort}.${pheno}*.png"), path("${cohort}.${pheno}*.csv")
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${exwas_regions_plot_script} \
          -c ${cohort} \
          -p ${pheno} \
          -g ${gene_file} \
          -s ${regions_sumstats} \
          -t ${pheno_table} \
          -m ${params.grouptest_maf} \
          -a "${params.grouptest_annotation}"
        """
    stub:
        """
        touch ${cohort}.${pheno}.regions.manhattan.png
        touch ${cohort}.${pheno}.regions.qq.png
        touch ${cohort}.${pheno}.exwas_regions.plots_manifest.csv
        """
}

process make_saige_gene_phewas_regions_plots {
    publishDir "${launchDir}/Plots/"
    label 'safe_to_skip', 'high_memory_plots'
    input:
        tuple val(cohort), val(chr), path(sumstats_file)
        path phenotype_descriptions
        path phewas_regions_plot_script
    output:
        tuple path('*.png'), path("${cohort}.${chr}.phewas_regions.plots_manifest.csv")
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${phewas_regions_plot_script} \
          -d ${phenotype_descriptions} \
          -s ${sumstats_file} \
          -c ${cohort} \
          -r ${chr}
        """
    stub:
        """
        touch phewas.png
        touch ${cohort}.${chr}.phewas_regions.plots_manifest.csv
        """
}

process make_exwas_report_methods_blurb {
    publishDir "${launchDir}/Summary"
    input:
        path exwas_methods_script
        val params_dict
    output:
        path('saige_exwas_methods.html')
    shell:
        """
        echo '${new JsonBuilder(params_dict).toPrettyString().replace(';', '|')}' > pipeline_params.txt
        ${params.my_python} ${exwas_methods_script}
        """
    stub:
        '''
        touch saige_exwas_methods.html
        '''
}

process make_exwas_report_src {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true
    input:
        path singles_plots
        path regions_plots
        path pheno_summary_plots
        path pheno_summary_table
        path singles_summary
        path regions_summary
        path methods_blurb
    output:
        path('src/', type: 'dir')
    shell:
        """
        mkdir src/

        for f in \$(echo "${singles_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        for f in \$(echo "${regions_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        for f in \$(echo "${pheno_summary_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
        echo \$f
        cp \$f src/
        done

        cp ${pheno_summary_table} src/
        cp ${singles_summary} src/
        cp ${regions_summary} src/

        cp ${methods_blurb} src/

        """
    stub:
        '''
        mkdir src
        touch src/test.png
        touch src/test.csv
        '''
}

process make_exwas_report {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true
    input:
        path report_source_dir
        path generate_html_script
    output:
        path('index.html')
        path('*.html')
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > regions_colnames.txt
        echo "${params.singles_col_names.collect().join('\n')}" > singles_colnames.txt
        ${params.my_python} ${generate_html_script} \
          -p ${params.bin_pheno_list.join(' ') + ' ' + params.quant_pheno_list.join(' ')} \
          -c ${params.cohort_list.join(' ')} \
          --pval ${params.p_cutoff_summarize} \
          --regionsColnames regions_colnames.txt \
          --singlesColnames singles_colnames.txt
        """
    stub:
        '''
        touch index.html
        '''
}

process collect_exwas_regions_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"
    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        for file_path in input_list:
            dfs.append(pd.read_csv(file_path))

        pd.concat(dfs,ignore_index=True).to_csv(output_file, index=False)
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}

process collect_exwas_singles_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"
    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()
                      if 'manifest.csv' in x.strip()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        with open(output_file, 'w') as f:
            f.write('filename\\n')
            for file_path in input_list:
                f.write(file_path + '\\n')
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}

process collect_gwas_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"
    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()
                      if 'manifest.csv' in x.strip()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        with open(output_file, 'w') as f:
            f.write('filename\\n')
            for file_path in input_list:
                f.write(file_path + '\\n')
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}

process make_gwas_report_src {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true
    input:
        path gwas_plots
        path pheno_summary_plots
        path pheno_summary_table
        path gwas_summary
    output:
        path('src/', type: 'dir')
    shell:
        """
        mkdir src/

        for f in \$(echo "${gwas_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
          cp \$f src/
        done

        for f in \$(echo "${pheno_summary_plots}" | sed 's|,||g' | sed 's|\\[||g' | sed 's|\\]||g')
        do
          cp \$f src/
        done

        cp ${pheno_summary_table} src/
        cp ${gwas_summary} src/
        """
    stub:
        '''
        mkdir src
        touch src/test.png
        touch src/test.csv
        '''
}

process make_gwas_report {
    publishDir "${launchDir}/Report/", mode: 'copy', overwrite: true
    input:
        path report_source_dir
        path generate_html_script
    output:
        path('index.html')
        path('*.html')
    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > gwas_colnames.txt
        ${params.my_python} ${generate_html_script} \
          -p ${params.bin_pheno_list.join(' ') + ' ' + params.quant_pheno_list.join(' ')} \
          -c ${params.cohort_list.join(' ')} \
          --pval ${params.p_cutoff_summarize} \
          --gwasColnames gwas_colnames.txt
        """
    stub:
        '''
        touch index.html
        '''
}

process collect_gene_phewas_regions_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"
    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        for file_path in input_list:
            dfs.append(pd.read_csv(file_path))

        pd.concat(dfs,ignore_index=True).to_csv(output_file, index=False)
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}

process make_saige_exwas_report {
    publishDir "${launchDir}", mode: 'copy'

    input:
        path all_plots, stageAs: 'Plots/*'
        path singles_summary
        path regions_summary
        path pheno_summaries
        path params_json
        path manifest_script
        path report_script
    output:
        path('SAIGE_ExWAS_Report.zip')
    shell:
        """
        ${params.my_python} ${manifest_script} \
            --params_json ${params_json} \
            --pheno_summaries ${pheno_summaries} \
            --singles_summary ${singles_summary} \
            --regions_summary ${regions_summary} \
            --plots_dir Plots/

        ${params.my_python} ${report_script} \
            --manifest saige_exwas_manifest.json \
            --output_zip SAIGE_ExWAS_Report.zip
        """
    stub:
        '''
        touch SAIGE_ExWAS_Report.zip
        '''
}

process make_saige_gene_phewas_report {
    publishDir "${launchDir}", mode: 'copy'

    input:
        path all_plots, stageAs: 'Plots/*'
        path singles_summary
        path regions_summary
        path pheno_summaries
        path params_json
        path plots_manifests, stageAs: 'PlotManifests/*'
        path manifest_script
        path report_script
    output:
        path('SAIGE_GenePheWAS_Report.zip')
    shell:
        """
        ${params.my_python} ${manifest_script} \
            --params_json ${params_json} \
            --pheno_summaries ${pheno_summaries} \
            --singles_summary ${singles_summary} \
            --regions_summary ${regions_summary} \
            --plots_dir Plots/ \
            --plots_manifests_dir PlotManifests/

        ${params.my_python} ${report_script} \
            --manifest saige_gene_phewas_manifest.json \
            --output_zip SAIGE_GenePheWAS_Report.zip
        """
    stub:
        '''
        touch SAIGE_GenePheWAS_Report.zip
        '''
}

process make_saige_variant_phewas_report {
    publishDir "${launchDir}", mode: 'copy'

    input:
        path all_plots, stageAs: 'Plots/*'
        path singles_summary
        path pheno_summaries
        path params_json
        path manifest_script
        path report_script
    output:
        path('SAIGE_VarPheWAS_Report.zip')
    shell:
        """
        ${params.my_python} ${manifest_script} \
            --params_json ${params_json} \
            --pheno_summaries ${pheno_summaries} \
            --singles_summary ${singles_summary} \
            --plots_dir Plots/

        ${params.my_python} ${report_script} \
            --manifest saige_variant_phewas_manifest.json \
            --output_zip SAIGE_VarPheWAS_Report.zip
        """
    stub:
        '''
        touch SAIGE_VarPheWAS_Report.zip
        '''
}

process make_saige_gwas_report {
    publishDir "${launchDir}", mode: 'copy'

    input:
        path all_plots, stageAs: 'Plots/*'
        path gwas_summary
        path pheno_summaries
        path params_json
        path manifest_script
        path report_script
    output:
        path('SAIGE_GWAS_Report.zip')
    shell:
        """
        ${params.my_python} ${manifest_script} \
            --params_json ${params_json} \
            --pheno_summaries ${pheno_summaries} \
            --gwas_summary ${gwas_summary} \
            --plots_dir Plots/

        ${params.my_python} ${report_script} \
            --manifest saige_gwas_manifest.json \
            --output_zip SAIGE_GWAS_Report.zip
        """
    stub:
        '''
        touch SAIGE_GWAS_Report.zip
        '''
}
