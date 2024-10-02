
process set_up_cohort {
    publishDir "${launchDir}/${cohort}/"
    machineType 'n2-standard-4'
    input:
        val cohort
        path cohort_script
        path pheno_covar_table
        path cohort_table
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(step2_fam, stageAs: 'Exome/*')
    output:
        tuple val(cohort), path('sample_list.txt')
        tuple val(cohort), path('saige_pheno_covars.txt')
    shell:
        """
        ${params.my_python} ${cohort_script} \
          --data ${pheno_covar_table} \
          --cohort ${cohort} \
          --samples ${cohort_table} \
          --step1Fam ${step1_fam} \
          --step2Fam ${step2_fam} \
          --id ${params.id_col}
        """
    stub:
        '''
        touch sample_list.txt
        touch saige_pheno_covars.txt
        '''
}

/*
process set_up_cohort_gwas_plink {
    publishDir "${launchDir}/${cohort}/"
    machineType 'n2-standard-4'
    input:
        val cohort
        path cohort_script
        path pheno_covar_table
        path cohort_table
        path(step1_fam, stageAs: 'Plink_QC/*')
        path(step2_fam, stageAs: 'Step2Fam/*')
    output:
        tuple val(cohort), path('sample_list.txt')
        tuple val(cohort), path('saige_pheno_covars.txt')
    shell:
        """
        ${params.my_python} ${cohort_script} \
          --data ${pheno_covar_table} \
          --cohort ${cohort} \
          --samples ${cohort_table} \
          --step1Fam ${step1_fam} \
          --step2Fam ${step2_fam} \
          --id ${params.id_col}
        """
    stub:
        '''
        touch sample_list.txt
        touch saige_pheno_covars.txt
        '''
}
*/

process set_up_cohort_gwas_bgen {
    publishDir "${launchDir}/${cohort}/"
    machineType 'n2-standard-4'
    input:
        val cohort
        path cohort_script
        path pheno_covar_table
        path cohort_table
        path(step1_fam, stageAs: 'Step1/*')
        path(bgen_samplefile, stageAs: 'bgenSample/*')
    output:
        tuple val(cohort), path('sample_list.txt')
        tuple val(cohort), path('saige_pheno_covars.txt')
    shell:
        """
        ${params.my_python} ${cohort_script} \
          --data ${pheno_covar_table} \
          --cohort ${cohort} \
          --samples ${cohort_table} \
          --step1Fam ${step1_fam} \
          --step2bgen_sample ${bgen_samplefile} \
          --id ${params.id_col}
        """
    stub:
        '''
        touch sample_list.txt
        touch saige_pheno_covars.txt
        '''
}

//TODO change exome to genetic data label for readability
process make_pheno_summaries {
    publishDir "${launchDir}/Summary/"
    machineType 'n2-standard-4'
    input:
        val cohort_list
        val bin_pheno_list
        val quant_pheno_list
        path(step1_fam, stageAs: 'Step1/*')
        path(step2_fam, stageAs: 'Step2/*')
        path pheno_covar_table
        path cohort_table
        path(pheno_table_script)
    output:
        path('pheno_summaries.csv')

    shell:
        """
        ${params.my_python} ${pheno_table_script} \
          -c ${cohort_list.join(' ')} \
          -b ${bin_pheno_list.join(' ')} \
          -q ${quant_pheno_list.join(' ')} \
          --step1Fam ${step1_fam} \
          --step2Fam ${step2_fam} \
          --data ${pheno_covar_table} \
          --samples ${cohort_table} \
          --id ${params.id_col}
        """
}
/*
process gpu_memory_reqs{
    machineType 'n2-standard-4'
    input:
        val(cohort)
        val(pheno)
        val()
    output:
    shell:
    '''

    '''

}
*/

/*
process check_bin_cohort_pheno_combo {
    machineType 'n2-standard-4'
    input:
        tuple val(cohort), val(pheno)
        val(pheno_table)
    output:
        tuple val(cohort), val(pheno), val(num_cases)
    script:
        // println ("${cohort}_${pheno}")
        pheno_count_file = new File(pheno_table.toString())
        rows = pheno_count_file.readLines().tail()*.split(',')
        rows = rows.findAll { it[0] == cohort }
        rows = rows.findAll { it[1] == pheno }
        match_row = rows.get(0)
        str_num_cases = match_row[4]
        // str_num_cases = rows.get(0)[4]
        // num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        """
        echo "${cohort} ${pheno}"
        """
}

process check_quant_cohort_pheno_combo {
    machineType 'n2-standard-4'
    input:
        tuple val(cohort), val(pheno)
        val(pheno_table)
    output:
        tuple val(cohort), val(pheno), val(num_cases)
    script:
        pheno_count_file = new File(pheno_table.toString())
        rows = pheno_count_file.readLines().tail()*.split(',')
        rows = rows.findAll { it[0] == cohort }
        rows = rows.findAll { it[1] == pheno }
        str_num_cases = rows.get(0)[2]
        num_cases = str_num_cases == '' ? 0 : str_num_cases.toDouble()
        '''
        '''
}
*/

// Actual process that parses the rows of the table:
process parse_pheno_summary_table {
    executor 'local'
    cache false
    input:
        path pheno_table_file // this is the appropriate path object
        val pheno_table_string // for executing Groovy code, must be val not path
    output:
        val list_of_info_lists
    script:
        exist_checks = 0
        sleep_interval = 100 // milliseconds
        max_sleep = 10 * 1000 // 10 seconds

        current_sleep = 0
        while (exist_checks < 3 && current_sleep <= max_sleep) {
            if (file(pheno_table_string).exists()) {
                exist_checks += 1
            }
            sleep(sleep_interval)
            current_sleep += sleep_interval
        }

        not_empty_checks = 0
        current_sleep = 0
        while (not_empty_checks < 3 && current_sleep <= max_sleep) {
            if ((new File(pheno_table_string)).empty()) {
                not_empty_checks += 1
            }
            sleep(sleep_interval)
            current_sleep += sleep_interval
        }

        assert file(pheno_table_string).exists() && !(new File(pheno_table_string)).empty() : 'File System Latency Issues. Please Try Re-Running'

        pheno_table_lines = (new File(pheno_table_string)).readLines()
        info_row_lists = []
        pheno_table_lines.each {
            line ->
            line_parts = line.toString().trim().replace('[', '').replace(']', '').split(',') as List
            info_row_lists.add(line_parts)
        }
        list_of_info_lists = [info_row_lists]

        '''
        echo "done"
        '''
}

include {
    paramToList
    get_script_file_names
} from '../processes/saige_helpers.nf'

include {
    make_pheno_covar_summary_plots
} from '../processes/saige_visualization.nf'

workflow SAIGE_PREPROCESSING {
    take:
        // all four of these are paths to input files
        pheno_covar_table
        cohort_table
        step1_fam
        step2_fam
        IS_PHEWAS
    main:
        MIN_BIN_CASES = params.min_bin_cases == null ? 50 : params.min_bin_cases
        MIN_QUANT_N = params.min_bin_cases == null ? 500 : params.min_quant_n

        cohort = Channel.fromList(params.cohort_list)
        bin_pheno_list = paramToList(params.bin_pheno_list)
        bin_pheno = Channel.fromList(bin_pheno_list)
        quant_pheno_list = paramToList(params.quant_pheno_list)
        quant_pheno = Channel.fromList(quant_pheno_list)
        chromosome = Channel.fromList(params.chromosome_list)

        script_name_dict = get_script_file_names()
        cohort_setup_script = script_name_dict['cohort_setup']

        if (params.ftype != null && params.ftype == 'BGEN') {
            (cohort_sample_lists, cohort_pheno_tables) = set_up_cohort_gwas_bgen(
            cohort, cohort_setup_script, pheno_covar_table,
            cohort_table, step1_fam, step2_fam)
        } else {
            (cohort_sample_lists, cohort_pheno_tables) = set_up_cohort(
            cohort, cohort_setup_script, pheno_covar_table,
            cohort_table, step1_fam, step2_fam)
        }

        pheno_table_script = script_name_dict['pheno_table']
        pheno_covar_plots_script = script_name_dict['pheno_covar_plots']

        // make pheno summary table, conditionally handle empty phenotype lists
        pheno_table = make_pheno_summaries(
                cohort.collect(),
                (bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
                (quant_pheno_list.size() == 0) ? '[]' : quant_pheno.toSortedList(),
                step1_fam, step2_fam,
                pheno_covar_table, cohort_table,
                pheno_table_script
                )

        if (!IS_PHEWAS) {
            pheno_plots = make_pheno_covar_summary_plots(
                cohort.collect(),
                (bin_pheno_list.size() == 0)  ? '[]' : bin_pheno.toSortedList(),
                (quant_pheno_list.size() == 0) ? '[]' : quant_pheno.toSortedList(),
                step1_fam, step2_fam,
                pheno_covar_table, cohort_table,
                pheno_covar_plots_script
                )
        }

        // parse the pheno-cohort summary table
        pheno_table_filename = "${launchDir}/Summary/pheno_summaries.csv"
        // parse the pheno-cohort summary table
        all_pheno_table_lines = parse_pheno_summary_table(pheno_table, pheno_table_filename)
        num_combos = params.cohort_list.size() * (bin_pheno_list.size() + quant_pheno_list.size())
        num_channel = Channel.fromList(1..num_combos)
        lines_with_num = num_channel.combine(all_pheno_table_lines)
        pheno_table_info = lines_with_num.map { i, all_lines ->
            all_lines.get(i) // Use i not i-1 because 0 is the header row
        }

        bin_pheno_info = pheno_table_info.filter { row -> bin_pheno_list.contains(row.get(1)) }
        quant_pheno_info = pheno_table_info.filter { row -> quant_pheno_list.contains(row.get(1)) }

        // sex-specific pheno handling
        sex_pheno_list = paramToList(params.sex_specific_pheno_file)

        // Filter out combinations with too few cases
        // For binary, must have MIN_BIN_CASES cases
        cohort_bin_pheno_case_ct = bin_pheno_info.map {
            row -> new Tuple(row.get(0), row.get(1), row.get(4) == '' ? 0 : row.get(4).toDouble())
        }

        keep_cohort_bin_pheno_combos = cohort_bin_pheno_case_ct \
            .filter {
            cohort, pheno, cases -> cases >= MIN_BIN_CASES } \
            .map { cohort, pheno, cases -> new Tuple(cohort, pheno) } \
            .filter {
                // sex-specific pheno handling
                // if the pheno is not in sex-specific list, keep it
                // OR if it IS in the list AND the cohort is in sex-stratified list, keep it
                cohort, pheno -> !sex_pheno_list.contains(pheno) || \
                (sex_pheno_list.contains(pheno) && params.sex_strat_cohort_list.contains(cohort))
            }

        // For quantitative, must have MIN_QUANT_N samples
        cohort_quant_pheno_ct = quant_pheno_info.map {
            row -> new Tuple(row.get(0), row.get(1), row.get(2) == '' ? 0 : row.get(2).toDouble())
        }

        keep_cohort_quant_pheno_combos = cohort_quant_pheno_ct \
            .filter { cohort, pheno, count -> count >= MIN_QUANT_N } \
            .map { cohort, pheno, count -> new Tuple(cohort, pheno) } \
            .filter {
                // sex-specific pheno handling
                // if the pheno is not in sex-specific list, keep it
                // OR if it IS in the list AND the cohort is in sex-stratified list, keep it
                cohort, pheno -> !sex_pheno_list.contains(pheno) || \
                (sex_pheno_list.contains(pheno) && params.sex_strat_cohort_list.contains(cohort))
            }

    emit:
        keep_cohort_bin_pheno_combos
        keep_cohort_quant_pheno_combos
        pheno_table
        cohort_sample_lists
        cohort_pheno_tables
}
//numGPU
