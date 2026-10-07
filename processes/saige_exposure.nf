/*
Processes for the phenotype exposure mode of the variant PheWAS workflow.
Each exposure is encoded as a dosage "marker" in a generated BGEN so that SAIGE Step 2
tests it against the Step 1 null models (which still model relatedness with the GRM).
*/

process make_exposure_bgen {
    publishDir "${launchDir}/Exposure_Genotypes/"
    input:
        path pheno_covar_table
        path cohort_table
        path(step1_fam, stageAs: 'Step1/*')
        path exposure_script
    output:
        tuple val('1'), path('exposures.chr1.{bgen,bgen.bgi}')
        path 'exposures.sample'
        path 'exposure_encoding.tsv'
        path 'exposure_cohort_sets.csv'
        path 'exposure_bgen.log'
    script:
        def bin_exposures = paramToList(params.bin_exposure_list)
        def cont_exposures = paramToList(params.cont_exposure_list)
        """
        ${params.my_python} ${exposure_script} \
          --data ${pheno_covar_table} \
          --id ${params.id_col} \
          --samples ${cohort_table} \
          --cohorts ${params.cohort_list.join(' ')} \
          --step1Fam ${step1_fam} \
          ${bin_exposures.size() > 0 ? '--bin_exposures ' + bin_exposures.join(' ') : ''} \
          ${cont_exposures.size() > 0 ? '--cont_exposures ' + cont_exposures.join(' ') : ''} \
          --cont_transform ${params.cont_exposure_transform} > exposure_bgen.log

        plink2 --vcf exposures.vcf dosage=DS \
          --export bgen-1.2 bits=8 ref-first \
          --out exposures.chr1 >> exposure_bgen.log

        ${params.my_bgenix} -g exposures.chr1.bgen -index -clobber
        """
    stub:
        """
        touch exposures.chr1.bgen
        touch exposures.chr1.bgen.bgi
        touch exposures.sample
        touch exposure_encoding.tsv
        touch exposure_cohort_sets.csv
        touch exposure_bgen.log
        """
}

process merge_and_filter_saige_exposure_phewas_output {
    publishDir "${launchDir}/${cohort_dir}/Sumstats/"
    // needs dynamic memory {} allocation
    maxRetries 5
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' } // Retry on OOM-related exit codes
    memory {
        def base_mem = params.host == 'AOU' ? 63.GB : 24.GB
        def attempt_mem = base_mem * task.attempt
        return attempt_mem
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), val(chr), path(chr_inputs), val(exposure)

        path merge_regions_script
        path exposure_encoding
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
          --exposure ${exposure} \
          --exposure_encoding ${exposure_encoding} \
          -s ${chr_inputs.join(' ')}
        """
    stub:
        """
        touch ${cohort_dir}.chr${chr}.variant_phewas.saige.gz
        touch ${cohort_dir}.chr${chr}.variant_phewas.filtered.saige.csv
        """
}

include {
    paramToList
} from '../processes/saige_helpers.nf'
