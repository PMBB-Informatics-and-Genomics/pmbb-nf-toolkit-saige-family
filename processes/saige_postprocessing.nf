
process merge_and_filter_saige_gene_regions_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)
        path merge_regions_script
        path pheno_summary_table
    output:
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_regions.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_regions.filtered.saige.csv")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_cauchy_tests.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_cauchy_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')} \
          --counts ${pheno_summary_table} \
          --regions
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.exwas_regions.saige.gz
        touch ${cohort_dir}.${pheno}.exwas_regions.filtered.saige.csv
        touch ${cohort_dir}.${pheno}.exwas_cauchy_tests.saige.gz
        touch ${cohort_dir}.${pheno}.exwas_cauchy_tests.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_singles_phewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(chr), path(chr_inputs)
        path merge_singles_script
    output:
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_singles.saige.gz")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_singles.filtered.saige.csv")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_ultrarare_tests.saige.gz")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_ultrarare_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_singles_script} \
          --chr ${chr} \
          --phewas \
          --rareVars \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        touch ${cohort_dir}.chr${chr}.gene_phewas_singles.saige.gz
        touch ${cohort_dir}.chr${chr}.gene_phewas_singles.filtered.saige.csv
        touch ${cohort_dir}.chr${chr}.gene_phewas_ultrarare_tests.saige.gz
        touch ${cohort_dir}.chr${chr}.gene_phewas_ultrarare_tests.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_regions_phewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(chr), path(pheno_inputs)
        path merge_regions_script
        path pheno_summary_table
    output:
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_regions.saige.gz")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_regions.filtered.saige.csv")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_cauchy_tests.saige.gz")
        tuple val(cohort_dir), val(chr), path("${cohort_dir}.chr${chr}.gene_phewas_cauchy_tests.filtered.saige.csv")
    shell:
    shell:
        """
        echo "${params.regions_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          --phewas \
          --regions \
          --chr ${chr} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          --counts ${pheno_summary_table} \
          -s ${pheno_inputs.join(' ')}
        """
    stub:
        """
        touch ${cohort_dir}.chr${chr}.gene_phewas_regions.saige.gz
        touch ${cohort_dir}.chr${chr}.gene_phewas_regions.filtered.saige.csv
        touch ${cohort_dir}.chr${chr}.gene_phewas_cauchy_tests.saige.gz
        touch ${cohort_dir}.chr${chr}.gene_phewas_cauchy_tests.filtered.saige.csv
        """
}

process merge_and_filter_saige_gene_singles_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)
        path merge_singles_script
    output:
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_singles.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_singles.filtered.saige.csv")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_ultrarare_tests.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.exwas_ultrarare_tests.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_singles_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')} \
          --rareVars
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.exwas_singles.saige.gz
        touch ${cohort_dir}.${pheno}.exwas_singles.filtered.saige.csv
        touch ${cohort_dir}.${pheno}.exwas_ultrarare_tests.saige.gz
        touch ${cohort_dir}.${pheno}.exwas_ultrarare_tests.filtered.saige.csv
        """
}

process merge_and_filter_saige_gwas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-16'
    memory '35GB'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr_list), path(chr_inputs)
        path merge_regions_script

    output:
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.gwas.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.${pheno}.gwas.filtered.saige.csv")

    shell:
        """
        echo "${params.gwas_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -p ${pheno} \
          -c colnames.txt \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}

        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.gwas.saige.gz
        touch ${cohort_dir}.${pheno}.gwas.filtered.saige.csv
        """
}

process gzipFiles {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-16'

    input:
    path(file_path)

    output:
    path "${file_path}.gz"
    script:
    """
    if [ -L "$file_path" ]; then
        actual_path=\$(readlink -f "$file_path")
        gzip "\$actual_path" > "${file_path}.gz"
        rm \$actual_path
    else
        gzip "$file_path" > "${file_path}.gz"
    fi

    """

    stub:
    """
    touch ${file_path}.gz
    """
}

process merge_and_filter_saige_variantphewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr), path(chr_inputs)

        path merge_regions_script
    output:
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.chr${chr}.variant_phewas.saige.gz")
        tuple val(cohort_dir), val(pheno), path("${cohort_dir}.chr${chr}.variant_phewas.filtered.saige.csv")
    shell:
        """
        echo "${params.singles_col_names.collect().join('\n')}" > colnames.txt
        cat colnames.txt
        ${params.my_python} ${merge_regions_script} \
          -c colnames.txt \
          --chr ${chr} \
          --phewas \
          --cohort ${cohort_dir} \
          --pvalue ${params.p_cutoff_summarize} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        touch ${cohort_dir}.chr${chr}.variant_phewas.saige.gz
        touch ${cohort_dir}.chr${chr}.variant_phewas.filtered.saige.csv
        """
}

process make_summary_regions_output {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'
    input:
        path(filtered_regions)
    output:
        path('saige_exwas_suggestive_regions.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${filtered_regions}".replace('[', '').replace(']', '').split()]
        output = 'saige_exwas_suggestive_regions.csv'

        for f in input_list:
            dfs.append(pd.read_csv(f))

        region_col = '${params.regions_col_names.keySet().toList().contains('Region') ? params.regions_col_names['Region'] : 'Region' }'
        pd.concat(dfs).sort_values(by=[region_col]).to_csv(output, index=False)
        """
    stub:
        '''
        touch saige_exwas_suggestive_regions.csv
        '''
}

process make_summary_singles_output {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'
    input:
        //stageAs: '?/*'
        path(filtered_singles)
    output:
        path('saige_exwas_suggestive_singles.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${filtered_singles}".replace("[", "").replace("]", "").split()]
        output = "saige_exwas_suggestive_singles.csv"

        for f in input_list:
            dfs.append(pd.read_csv(f))

        chr_col = '${params.singles_col_names.keySet().toList().contains('CHR') ? params.singles_col_names['CHR'] : 'CHR' }'
        pos_col = '${params.singles_col_names.keySet().toList().contains('POS') ? params.singles_col_names['POS'] : 'POS' }'
        pd.concat(dfs).sort_values(by=[chr_col, pos_col]).to_csv(output, index=False)
        """
    stub:
        '''
        touch saige_exwas_suggestive_singles.csv
        '''
}

process make_summary_suggestive_gwas {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'
    input:
        path(filtered_singles)
    output:
        path('saige_gwas_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${filtered_singles}".replace("[", "").replace("]", "").split()]
        output = "saige_gwas_suggestive.csv"

        for f in input_list:
            dfs.append(pd.read_csv(f))

        chr_col = '${params.gwas_col_names.keySet().toList().contains('CHR') ? params.gwas_col_names['CHR'] : 'CHR' }'
        pos_col = '${params.gwas_col_names.keySet().toList().contains('POS') ? params.gwas_col_names['POS'] : 'POS' }'
        pd.concat(dfs).sort_values(by=[chr_col, pos_col]).to_csv(output, index=False)
        """
    stub:
        '''
        touch saige_gwas_suggestive.csv
        '''
}

//, arity: '1..*'
process gwas_make_biofilter_positions_input {
    publishDir "${launchDir}/Annotations/"
    machineType 'n2-standard-4'
    input:
        //, stageAs: "?/*"
        path(filtered_sumstats)
    output:
        path('gwas_biofilter_input_positions.txt')
    script:
        """
        #! ${params.my_python}
        import pandas as pd

        dfs = []
        input_list = '${filtered_sumstats.join(' ')}'.split()
        print(input_list)
        for f in input_list:
            temp = pd.read_table(f, sep=",")
            print(temp.head())
            temp['chromosome'] = temp['chromosome'].astype(str).str.replace('chr', '').astype(int)
            temp = temp.rename(columns={'chromosome': 'CHR', 'base_pair_location': 'POS', 'p_value': 'P', 'variant_id': 'ID'})
            dfs.append(temp)
        all = pd.concat(dfs)
        keep_cols = ['CHR', 'ID', 'POS']
        all[keep_cols].to_csv('gwas_biofilter_input_positions.txt', header=False, index=False, sep=' ')
        """
    stub:
        '''
        touch gwas_biofilter_input_positions.txt
        '''
}

process make_summary_table_with_annot {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'
    input:
        path(all_filtered_sumstats)
        tuple val(data_nickname), path(biofilter_annots)
    output:
        path('saige_gwas_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "saige_gwas_suggestive.csv"

        for f in input_list:
            dfs.append(pd.read_csv(f))

        all_meta = pd.concat(dfs).sort_values(by='variant_id')

        annot_df = pd.read_csv('${biofilter_annots}', index_col='Var_ID')
        all_meta[['Gene', 'RSID']] = annot_df.loc[all_meta['variant_id'], ['Gene', 'RSID']].values

        all_meta.to_csv(output, index=False)
        """
    stub:
        '''
        touch saige_gwas_suggestive.csv
        '''
}
