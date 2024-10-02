process filter_chr_group_file {
    publishDir "${launchDir}/PheWAS_GroupFiles/"
    machineType 'n2-standard-4'
    input:
        tuple val(chr), path(chr_group_file)
        path(gene_list_file)
    output:
        tuple val(chr), path("filtered_group_file.${chr}.txt")
    script:
        """
        #! ${params.my_python}

        keep_genes = open('${gene_list_file}').read().splitlines()

        keep_lines = []

        with open('${chr_group_file}') as f:
            for line in f:
                if line.split('\\t', 2)[0] in keep_genes:
                    keep_lines.append(line.strip('\\n'))

        open('filtered_group_file.${chr}.txt', 'w+').write('\\n'.join(keep_lines) + '\\n')
        """
    stub:
        """
        touch filtered_group_file.${chr}.txt
        """
}

process call_saige_gene_step2_bin {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.gz"), path("${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz"), path("${cohort_dir}.${pheno}.${chr}.markerList.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --is_Firth_beta=${params.use_firth} \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr} \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}
        gzip -9 ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.markerList.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        touch ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.markerList.txt.gz
        """
}

process call_saige_gene_step2_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.gz"), path("${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz"), path("${cohort_dir}.${pheno}.${chr}.markerList.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}
        gzip -9 ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.markerList.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        touch ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.markerList.txt.gz
        """
}

process call_saige_gene_step2_bin_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"
    machineType 'n2-standard-4'
    input:
        // variables, step 1 outputs
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.gz"), path("${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz"), path("${cohort_dir}.${pheno}.${chr}.markerList.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --is_Firth_beta=${params.use_firth} \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr} \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}
        gzip -9 ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.markerList.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        touch ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.markerList.txt.gz
        """
}

process call_saige_gene_step2_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Gene_Results/"
    machineType 'n2-standard-4'
    input:
        // variables, step 1 outputs
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        path(chr_group_file)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.gz"), path("${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz"), path("${cohort_dir}.${pheno}.${chr}.markerList.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --maxMAF_in_groupTest="${params.grouptest_maf}" \
         --annotation_in_groupTest="${params.grouptest_annotation}" \
         --groupFile=${chr_group_file} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --is_output_markerList_in_groupTest=TRUE \
         --is_single_in_groupTest=TRUE \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}
        gzip -9 ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.markerList.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        touch ${cohort_dir}.${pheno}.${chr}.singleAssoc.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.markerList.txt.gz
        """
}

include {
    paramToList
} from '../processes/saige_helpers.nf'

workflow SAIGE_GENE_STEP2 {
    take:
        step1_bin_output
        step1_quant_output
        use_step2_prefix
        IS_PHEWAS
    main:
        cohort = Channel.fromList(params.cohort_list)
        bin_pheno = Channel.fromList(paramToList(params.bin_pheno_list))
        quant_pheno = Channel.fromList(paramToList(params.quant_pheno_list))
        chromosome = Channel.fromList(params.chromosome_list)

        plink_suffixes_list = ['.bed', '.bim', '.fam']

        /*
        Step 1 -> Step 2 Channel Emission Tuples
        Step 1 Out:  cohort, phenotype, rda, var
        Combine:     chromosome
        Step 2 In:   cohort, phenotype, rda, var, chromosome
        */
        step2_bin_input = step1_bin_output.combine(chromosome)
        step2_quant_input = step1_quant_output.combine(chromosome)

        // Lets make the paths for the chromosome-separated group files
        chr_group_files = chromosome.map { chr -> new Tuple(chr, params.group_file_prefix + chr + '.txt') }

        if (IS_PHEWAS) {
            // For a PheWAS, it helps speed things up a lot if we filter the group files
            gene_list = params.gene_list_file

            // Replace the chr_group_files channel with a channel that has the same values
            // But now with paths to the filtered files
            chr_group_files = filter_chr_group_file(chr_group_files, gene_list)
        }

        // For the join, we need to combine our group files with cohort and phenotype
        group_files_parallel_bin = chr_group_files.combine(cohort).combine(bin_pheno).map {
            chr, group_file, cohort, pheno -> new Tuple(cohort, pheno, chr, group_file)
        }
        group_files_parallel_quant = chr_group_files.combine(cohort).combine(quant_pheno).map {
            chr, group_file, cohort, pheno -> new Tuple(cohort, pheno, chr, group_file)
        }

        // This attaches our group files to our step2 input by joining on (cohort, pheno, chr)
        chr_group_file_bin = step2_bin_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(group_files_parallel_bin, by: [0, 1, 2]).map {
            cohort, pheno, chr, group_file -> group_file
        }

        chr_group_file_quant = step2_quant_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(group_files_parallel_quant, by: [0, 1, 2]).map {
            cohort, pheno, chr, group_file -> group_file
        }

        // Now, we make sure the plink files are staged for each job
        exome_plink_file_bin = step2_bin_input.map {
            cohort, pheno, rda, var, chr -> \
            new Tuple(*plink_suffixes_list.collect { ext -> use_step2_prefix + ext })
        }
        exome_plink_file_quant = step2_quant_input.map {
            cohort, pheno, rda, var, chr -> \
            new Tuple(*plink_suffixes_list.collect { ext -> use_step2_prefix + ext })
        }

        // Finally, we can call SAIGE Step 2!
        if (params.use_sparse_GRM) {
            // Define sparse GRM input
            sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)
            // Binary SAIGE-GENE step 2 with GRM
            (step2_bin_output, step2_bin_logs) = call_saige_gene_step2_bin_with_sparse_GRM(
                step2_bin_input, chr_group_file_bin,
                exome_plink_file_bin, sparse_grm_input)
            // Quantitative SAIGE-GENE step 2 with GRM
            (step2_quant_output, step2_quant_logs) = call_saige_gene_step2_quant_with_sparse_GRM(
                step2_quant_input, chr_group_file_quant,
                exome_plink_file_quant, sparse_grm_input)
        }
        else {
            // Binary SAIGE-GENE step 2 without pre-computed GRM
            (step2_bin_output, step2_bin_logs) = call_saige_gene_step2_bin(
                step2_bin_input, chr_group_file_bin, exome_plink_file_bin)
            // Quantitative SAIGE-GENE step 2 without pre-computed GRM
            (step2_quant_output, step2_quant_logs) = call_saige_gene_step2_quant(
                step2_quant_input, chr_group_file_quant, exome_plink_file_quant)
        }
    emit:
        step2_bin_output
        step2_quant_output
}
