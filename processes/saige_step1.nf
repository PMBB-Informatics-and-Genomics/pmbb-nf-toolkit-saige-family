process plink_qc_for_step1 {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'step1_plink_qc'
    cpus 16
    memory { params.host == 'AOU' ? '63GB' : '24GB' }

    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    script:
        // use_mem = params.host == 'AOU' ? '62000' : '23000'
        use_mem = Math.floor(task.memory.toMega() * 0.85) // 85% of allocated memory in MB
        is_max_vars_null = params.max_vars_for_GRM == null ? 'YES' : 'NO'
        use_max_vars = params.max_vars_for_GRM == null ? '150000' : params.max_vars_for_GRM
        use_min_vars = params.min_vars_for_GRM == null ? '30000' : params.min_vars_for_GRM
        use_prune_r2 = params.pruning_r2_for_GRM == null ? '0.6' : params.pruning_r2_for_GRM
        // max_vars_for_GRM = params.max_vars_for_GRM == null ? '150000' : params.max_vars_for_GRM
        // min_vars_for_GRM = params.min_vars_for_GRM == null ? '30000' : params.min_vars_for_GRM
        """
        num_samples=\$(< ${sample_list} wc -l)
        rare_count=\$(awk -v n="\$num_samples" -v maf="${params.maf}" 'BEGIN { printf "%d", n * 2 * maf }')

        plink2 --indep-pairwise 50kb 1 ${use_prune_r2} \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --geno ${params.geno} \
            --maf ${params.maf} \
            --not-chr X,Y,MT,XY,PAR1,PAR2 \
            --hwe ${params.hwe} \
            --snps-only \
            --out step1_pruning
        
        num_maf=\$(< step1_pruning.prune.in wc -l)

        echo "Found \${num_maf} independent SNPs after filtering (min: ${use_min_vars}, max: ${use_max_vars})"
        
        # if the number of variants is less than the minimum, exit
        if [ "\$num_maf" -lt ${use_min_vars} ]; then
            echo "ERROR: Insufficient variants for stable GRM (\${num_maf} < ${use_min_vars})"
            echo "Please reduce GRM MAF threshold or provide more markers"
            exit 1
        
        # if thenumber of variants is greater than the maximum, and no maximum 
        # is set, then force user to explicitly set a maximum rather than silently reducing number of variants
        elif [ "${is_max_vars_null}" = "YES" ] && [ "\$num_maf" -gt 150000 ]; then
            echo "Error: Greater than 150k variants left which could cause step 1 to run for a long time."
            echo "You have 2 options: (1) set params.max_vars_for_GRM to desired number of variants, and we will randomly subset them."
            echo "Or (2) reduce number of markers"
            exit 1
        
        # if the number of variants is less than the maximum, and no maximum is set, then use all variants
        elif [ "${is_max_vars_null}" = "YES" ] && [ "\$num_maf" -lt 150000 ]; then
            echo "INFO: \$num_maf < 150000. Using all \${num_maf} variants for GRM construction"
            cp step1_pruning.prune.in  step1.markerid.list
        
        # if the number of variants is less than or equal to the maximum, and a maximum is set, then use entire set
        elif [ "${is_max_vars_null}" = "NO" ] && [ "\$num_maf" -lt ${use_max_vars} ] || [ "\$num_maf" -eq ${use_max_vars} ]; then
            echo "INFO: \$num_maf <= ${use_max_vars}. Using all \${num_maf} variants for GRM construction"
            cp step1_pruning.prune.in  step1.markerid.list
        
        # if the number of variants is greater than the maximum, and a maximum is set, then randomly subset down to maximum
        elif [ "${is_max_vars_null}" = "NO" ] && [ "\$num_maf" -gt ${use_max_vars} ]; then
            echo "INFO: \$num_maf > ${use_max_vars}, randomly selecting ${use_max_vars} variants"
            shuf -n ${use_max_vars} step1_pruning.prune.in > step1.markerid.list
        fi

        # Copy the filtered IDs to the output file
        echo "Filtered IDs written to step1.markerid.list"

        sort step1.markerid.list | uniq > step1.markerid.dedup.txt

        stdbuf -e0 -o0 plink2 --make-bed \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --extract step1.markerid.dedup.txt \
            --out step1_plink > plink_qc.log
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

process plink_qc_for_step1_saige_gene {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'step1_plink_qc'
    cpus 16
    // needs dynamic memory {} allocation
    maxRetries 5 // Retry up to 5 times
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    memory {
        def base_mem = params.host == 'AOU' ? 63.GB : 24.GB
        def attempt_mem = base_mem * task.attempt
        return attempt_mem
    }
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    script:
        // use_mem = params.host == 'AOU' ? '62000' : '23000'
        use_mem = Math.floor(task.memory.toMega() * 0.85) // 95% of allocated memory in MB
        is_max_vars_null = params.max_vars_for_GRM == null ? 'YES' : 'NO'
        use_max_vars = params.max_vars_for_GRM == null ? '150000' : params.max_vars_for_GRM
        use_min_vars = params.min_vars_for_GRM == null ? '30000' : params.min_vars_for_GRM
        use_min_rare_vars = params.min_rare_vars_for_GRM == null ? '300' : params.min_rare_vars_for_GRM
        use_prune_r2 = params.pruning_r2_for_GRM == null ? '0.6' : params.pruning_r2_for_GRM
        """
        num_samples=\$(< ${sample_list} wc -l)
        rare_count=\$(awk -v n="\$num_samples" -v maf="${params.maf}" 'BEGIN { printf "%d", n * 2 * maf }')

        plink2 --indep-pairwise 50kb 1 ${use_prune_r2} \
            --freq counts \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --geno ${params.geno} \
            --not-chr X,Y,MT,XY,PAR1,PAR2 \
            --hwe ${params.hwe} \
            --snps-only \
            --out step1_pruning
        
        wc -l step1_pruning.acount

        awk '{ if ((\$5 < (\$6-\$5) && \$5 >= 10 && \$5 <= 20) || (\$5 > (\$6-\$5) && \$6-\$5 >= 10 && \$6-\$5 <= 20)) { print \$2} }' ./step1_pruning.acount > temp_mac10to20_ids.txt
        awk '{ if ((\$5 < (\$6-\$5) && \$5 > 20) || (\$5 > (\$6-\$5) && \$6-\$5 > 20)) { print \$2} }' step1_pruning.acount > temp_mac20ormore_ids.txt
        awk -v ct="\$rare_count" '{ if ((\$5 < (\$6-\$5) && \$5 >= ct) || (\$5 > (\$6-\$5) && \$6-\$5 >= ct)) { print \$2} }' step1_pruning.acount > temp_maf_min.txt

        comm -12 <(sort temp_maf_min.txt) <(sort step1_pruning.prune.in) > prune_and_maf.txt

        num_mac10=\$(< temp_mac10to20_ids.txt wc -l)
        num_mac20=\$(< temp_mac20ormore_ids.txt wc -l)
        num_maf=\$(< prune_and_maf.txt wc -l)

        echo "\${num_mac10} Mac 10-20 \${num_mac20} MAC >20 \${num_maf} variants > MAF."

        if [ "\$num_mac10" -lt ${use_min_rare_vars} ]; then
            echo "Error: Less than 300 SNPs with MAC between 10 and 20"
            echo "Please provide more rare markers or set params.use_min_rare_vars to a string less than ${use_min_rare_vars}"
            exit 1
        elif [ "\$num_mac20" -lt ${use_min_rare_vars} ]; then
            echo "Error: Less than 300 SNPs with MAC greater than 20"
            echo "This is highly unusual. Please check your data"
            exit 1
        
        elif [ "\$num_maf" -lt ${use_min_vars} ]; then
            echo "ERROR: Insufficient variants for stable GRM (\${num_maf} < ${use_min_vars})"
            echo "Please reduce GRM MAF threshold or provide more markers"
            exit 1
        
        # if the number of variants is greater than the maximum, and no maximum 
        # is set, then force user to explicitly set a maximum rather than silently reducing number of variants
        elif [ "${is_max_vars_null}" = "YES" ] && [ "\$num_maf" -gt 150000 ]; then
            echo "Error: Greater than 150k variants left which could cause step 1 to run for a long time."
            echo "You have 2 options: (1) set params.max_vars_for_GRM to desired number of variants, and we will randomly subset them."
            echo "Or (2) reduce number of markers"
            exit 1
        
        # if the number of variants is less than the maximum, and no maximum is set, then use all variants
        elif [ "${is_max_vars_null}" = "YES" ] && [ "\$num_maf" -lt 150000 ]; then
            echo "INFO: \$num_maf < 150000. Using all \${num_maf} variants for GRM construction"
            cp prune_and_maf.txt  step1.markerid.list
        
        # if the number of variants is less than or equal to the maximum, and a maximum is set, then use entire set
        elif [ "${is_max_vars_null}" = "NO" ] && [ "\$num_maf" -lt ${use_max_vars} ] || [ "\$num_maf" -eq ${use_max_vars} ]; then
            echo "INFO: \$num_maf <= ${use_max_vars}. Using all \${num_maf} variants for GRM construction"
            cp prune_and_maf.txt  step1.markerid.list
        
        # if the number of variants is greater than the maximum, and a maximum is set, then randomly subset down to maximum
        elif [ "${is_max_vars_null}" = "NO" ] && [ "\$num_maf" -gt ${use_max_vars} ]; then
            echo "INFO: \$num_maf > ${use_max_vars}, randomly selecting ${use_max_vars} variants"
            shuf -n ${use_max_vars} prune_and_maf.txt > step1.markerid.list
        fi

        # Copy the filtered IDs to the output file
        shuf -n 300 temp_mac10to20_ids.txt >> ./step1.markerid.list
        shuf -n 300 temp_mac20ormore_ids.txt >> ./step1.markerid.list
        echo "Filtered IDs written to step1.markerid.list"

        sort step1.markerid.list | uniq > step1.markerid.dedup.txt

        stdbuf -e0 -o0 plink2 --make-bed \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --extract step1.markerid.dedup.txt \
            --out step1_plink > plink_qc.log
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}


String get_covar_list_args(String cohort, cohort_cat_covars, cohort_cont_covars) {
    String output = ''

    // if there are either continuous or categorical covariates, then add them to the covarColList
    if (cohort_cont_covars.size() > 0 || cohort_cat_covars.size() > 0) {
        output += '--covarColList='
        
        // Collect all covariates in a list first
        def all_covars = []
        if (cohort_cont_covars.size() > 0) {
            all_covars.addAll(cohort_cont_covars)
        }
        if (cohort_cat_covars.size() > 0) {
            all_covars.addAll(cohort_cat_covars)
        }
        
        output += all_covars.join(',') + ' '
    }

    // if there are categorical covariates, then add them to the qCovarColList
    // note that qCovarColList is a list of categorical covariates, while covarColList is a list of all covariates (continuous and categorical)
    if (cohort_cat_covars.size() > 0) {
        output += '--qCovarColList=' + cohort_cat_covars.join(',') + ' '
    }

    return output
}

process dnanexus_skip_plink_qc_for_step1 {
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)

    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'

    script:
    """
    # Check if the required PLINK files exist in the launch directory
    if [[ ! -f "${plink_prefix}.bed" || ! -f "${plink_prefix}.bim" || ! -f "${plink_prefix}.fam" ]]; then
        echo "Error: One or more PLINK files (.bed, .bim, .fam) are missing in the launch directory." >&2
        exit 1
    fi

    # Create symbolic links to mimic output files
    ln -s "${plink_prefix}.bed" step1_plink.bed
    ln -s "${plink_prefix}.bim" step1_plink.bim
    ln -s "${plink_prefix}.fam" step1_plink.fam

    # Create a dummy log file
    echo "PLINK file check completed successfully." > plink_qc.log
    """

    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

process dnanexus_skip_step1_quant {
    //for this to work, output files must exist in the launch dir
    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    script:
    """
    # Check if the required PLINK files exist in the launch directory
    if [[ ! -f "${pheno}.varianceRatio.txt" || ! -f "${pheno}.rda" || ! -f "${pheno}.log" ]]; then
        echo "Error: One or more PLINK files (.bed, .bim, .fam) are missing in the launch directory." >&2
        exit 1
    fi

    # Create symbolic links to mimic output files
    ln -s "${pheno}.varianceRatio.txt" ${pheno}.varianceRatio.txt
    ln -s "${pheno}.rda" ${pheno}.rda
    ln -s "${pheno}.log" ${pheno}.log

    # Create a dummy log file
    echo "STEP 1 file check completed successfully." > step1.log
    """

    stub:
        """
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.rda
        touch ${pheno}.log
        touch step1.log
        """
}

process dnanexus_skip_step1_bin {
    input:
    tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)

    output:
    tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
    path "${pheno}.log"

    script:
    """
    # Use Nextflow's launchDir
    if [[ ! -f "${launchDir}/${pheno}.varianceRatio.txt" || ! -f "${launchDir}/${pheno}.rda" || ! -f "${launchDir}/${pheno}.log" ]]; then
        echo "Error: One or more files are missing in ${launchDir}" >&2
        ls -la "${launchDir}"
        exit 1
    fi

    ln -sf "${launchDir}/${pheno}.varianceRatio.txt" ./${pheno}.varianceRatio.txt
    ln -sf "${launchDir}/${pheno}.rda" ./${pheno}.rda
    ln -sf "${launchDir}/${pheno}.log" ./${pheno}.log

    echo "STEP 1 file check completed successfully." > step1.log
    """
}

//machineType { (params.GPU == 'ON') ? getGpuMachineTypeChannel(plink_set[1].readLines().size, pheno_file.readLines().size) : 'n2-standard-4' }

process call_saige_step1_bin {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    memory { params.host == 'AOU' ? '63GB' : '24GB' }
    cpus 16
    
    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        // use_mem = params.host == 'AOU' ? '56' : '16'
        use_mem = Math.floor(task.memory.toGiga() * 0.85) // 95% of allocated memory in GB
        // Calculate nThreads as task.cpus - 1, with a minimum of 1
        use_threads = Math.max(1, task.cpus - 1)
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        """
        if [ "${params.GPU}" = "ON" ] && [ "${params.host}" = "LPC" ]; then
            mpirun -n 8 singularity exec --bind /path/to/data/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "DNAnexus" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno}\
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
         elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "AOU" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno}\
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --memoryChunk=${use_mem} \
            --nThreads=${use_threads} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --LOCO=${params.LOCO} \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi

        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_quant {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    // errorStrategy 'retry'
    memory { params.host == 'AOU' ? '63GB' : '24GB' }

    cpus 16

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        // Calculate memory 
        use_mem = Math.floor(task.memory.toGiga() * 0.85) // 95% of allocated memory in GB
        // Calculate nThreads as task.cpus - 1, with a minimum of 1
        use_threads = Math.max(1, task.cpus - 1)
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        """
        if [ "${params.GPU}"  = "ON" ]; then
            mpirun -n 8 singularity exec --bind /path/to/data/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=quantitative \
            --outputPrefix=${pheno} \
            --invNormalize=TRUE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=quantitative \
            --outputPrefix=${pheno} \
            --memoryChunk=${use_mem} \
            --nThreads=${use_threads} \
            --LOCO=${params.LOCO} \
            --invNormalize=TRUE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi

        """
    stub:
        """
        touch ${pheno}.log
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        """
}

process call_saige_step1_survival {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    // errorStrategy 'retry'
    memory { params.host == 'AOU' ? '63GB' : '24GB' }

    cpus 16

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        use_mem = Math.floor(task.memory.toGiga() * 0.85) // 95% of allocated memory in GB
        // Calculate nThreads as task.cpus - 1, with a minimum of 1
        use_threads = Math.max(1, task.cpus - 1)
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        """
        if [ "${params.GPU}" = "ON" ] && [ "${params.host}" = "LPC" ]; then
            mpirun -n 8 singularity exec --bind /path/to/data/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --nThreads=${use_threads} \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "DNAnexus" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --memoryChunk=${use_mem} \
            --nThreads=14 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            ${vr_flag} \
            --minMAFforGRM=${params.maf} \
            --maxMissingRateforGRM=${params.geno} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --LOCO=${params.LOCO} \
            --nThreads=${use_threads} \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_bin_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    cpus 16

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        use_threads = Math.max(1, task.cpus - 1)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --useSparseGRMforVarRatio=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         ${vr_flag} \
         --sampleIDColinphenoFile=${params.id_col} \
         --traitType=binary \
         --outputPrefix=${pheno} \
         --nThreads=${use_threads} \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    cpus 16

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        use_threads = Math.max(1, task.cpus - 1)
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --useSparseGRMforVarRatio=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         ${vr_flag} \
         --sampleIDColinphenoFile=${params.id_col} \
         --outputPrefix=${pheno} \
         --traitType=quantitative \
         --invNormalize=TRUE \
         --nThreads=${use_threads} \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_survival_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/", enabled: {params.host != 'AOU'}
    label 'saige_step1'
    label 'saige_process'
    cpus 16

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
        val is_gene
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        use_threads = Math.max(1, task.cpus - 1)
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        vr_flag = is_gene ? '--isCateVarianceRatio=TRUE' : '--isCateVarianceRatio=FALSE'
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --useSparseGRMforVarRatio=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         ${vr_flag} \
         --sampleIDColinphenoFile=${params.id_col} \
         --outputPrefix=${pheno} \
         --traitType=survival \
         --eventTimeCol=${params.event_time_col} \
         --eventTimeBinSize=${params.event_time_bin} \
         --invNormalize=TRUE \
         --nThreads=${use_threads} \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

include {
    paramToList
    getGpuMachineTypeChannel
} from '../processes/saige_helpers.nf'

workflow SAIGE_STEP1 {
    take:
        cohort_sample_lists
        cohort_pheno_tables
        keep_cohort_bin_pheno_combos
        keep_cohort_quant_pheno_combos
        keep_cohort_survival_pheno_combos
        use_step1_prefix
        IS_GENE
    main:
        bin_pheno_list = paramToList(params.bin_pheno_list)
        bin_pheno = Channel.fromList(bin_pheno_list)
        quant_pheno_list = paramToList(params.quant_pheno_list)
        quant_pheno = Channel.fromList(quant_pheno_list)
        survival_pheno_list = paramToList(params.survival_pheno_list)
        survival_pheno = Channel.fromList(survival_pheno_list)
        plink_suffixes_list = ['.bed', '.bim', '.fam']

        // First decide how to set up plink files - no GRM means we have to run QC. with GRM, we can use them directly
        if (params.use_sparse_GRM) {
            // We can use the input plink files directly - they should be the set of variants used to make the GRM
            cohort_plinkset = cohort_sample_lists.map { cohort, sample_list -> \
                new Tuple(cohort, new Tuple(*plink_suffixes_list.collect { ext -> use_step1_prefix + ext }))
            }
        }
        else {
            // set up inputs for plink_qc
            plink_inputs_tuple = new Tuple(*plink_suffixes_list.collect { ext -> use_step1_prefix + ext })

            if (IS_GENE) {
                // call plink_qc for each cohort (gene version, samples rare variants)
                //(cohort_plinkset, logs) = plink_qc_for_step1_saige_gene(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
                (cohort_plinkset, logs) = plink_qc_for_step1_saige_gene(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
            }
            else {
                if (params.host == 'DNAnexus') {
                    // Check for existing PLINK files instead of running QC
                    (cohort_plinkset, logs) = dnanexus_skip_plink_qc_for_step1(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
                }
                else {
                    // call plink_qc for each cohort
                    (cohort_plinkset, logs) = plink_qc_for_step1(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
                }
            }
        }

        // Add the sample lists and phenotype tables in so we get:
        // (cohort, (plink_trio), sample_list, pheno_table)
        cohort_plinkset_sample_pheno = cohort_plinkset.join(cohort_sample_lists).join(cohort_pheno_tables)

        /*
        Step 1 Input Emission Tuples
        Need: (cohort, phenotype, (plink set), sample_file, pheno_file)
        Have:
        (cohort, (plink_set), sample_file, pheno_file)
        (phenotype)
        Steps:
        1. .combine(bin_pheno or quant_pheno) -> (cohort, (plink set), sample_file, pheno_file, phenotype)
        2. .map() -> rearrange to get cohort, pheno first
        3. .join() -> keep only cohort/pheno combos with enough sample size
        */
        step1_bin_input = cohort_plinkset_sample_pheno \
            .combine(bin_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_bin_pheno_combos, by: [0, 1])

        step1_quant_input = cohort_plinkset_sample_pheno \
            .combine(quant_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_quant_pheno_combos, by: [0, 1])

        step1_survival_input = cohort_plinkset_sample_pheno \
            .combine(survival_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_survival_pheno_combos, by: [0, 1])

        // Now figure out which version of Step 1 We're Calling
        if (params.use_sparse_GRM) {
            sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)

            (step1_bin_output, logs) = call_saige_step1_bin_with_sparse_GRM(step1_bin_input, sparse_grm_input, IS_GENE)
            (step1_quant_output, logs) = call_saige_step1_quant_with_sparse_GRM(step1_quant_input, sparse_grm_input, IS_GENE)
            (step1_surival_output, logs) = call_saige_step1_survival_with_sparse_GRM(step1_survival_input, sparse_grm_input, IS_GENE)
        } else {
            if (params.host == 'DNAnexus') {
                (step1_bin_output, logs) = dnanexus_skip_step1_bin(step1_bin_input)
            }
            //lpc or aou run
            else {
                (step1_bin_output, logs) = call_saige_step1_bin(step1_bin_input, IS_GENE)
            }
            //if (params.host == "DNAnexus") {
            //    (step1_bin_output, logs) = dnanexus_skip_step1_quant(step1_quant_input)
            //}
            //lpc or aou run
            //else{
            (step1_quant_output, logs) = call_saige_step1_quant(step1_quant_input, IS_GENE)
            //}
            (step1_surival_output, logs) = call_saige_step1_survival(step1_survival_input, IS_GENE)
        }
    emit:
        step1_bin_output
        step1_quant_output
        step1_surival_output
}
