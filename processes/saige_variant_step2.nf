/*
This is a collection of SAIGE Step 2 processes
*/

process filter_snps_plink {
    publishDir "${launchDir}/Filtered_Files/"
    machineType 'n2-standard-4'
    input:
        tuple val(chr), path(plink_set)
        path(snp_list_file)
    output:
        tuple val(chr), path("phewas_input.chr${chr}.{bed,bim,fam}")
    shell:
        """
        plink --bfile ${plink_set[0].toString().split('/')[-1].replace('.bed', '')} \
        --extract range ${snp_list_file} \
        --make-bed \
        --out phewas_input.chr${chr}
        """
    stub:
        """
        touch phewas_input.chr${chr}.bed
        touch phewas_input.chr${chr}.bim
        touch phewas_input.chr${chr}.fam
        """
}

process filter_snps_bgen {
    publishDir "${launchDir}/Filtered_Files/"
    machineType 'n2-standard-4'
    input:
        tuple val(chr), path(bgen_set)
        path(snp_list_file)
    output:
        tuple val(chr), path("phewas_input.chr${chr}.{bgen,bgen.bgi}")
    shell:
        """
        awk '{print \$1 ":" \$2 "-" \$3}' ${snp_list_file} > bgen_formatted_snplist.txt
        ${params.my_bgenix} -g ${bgen_set[0].toString()} -incl-range bgen_formatted_snplist.txt > phewas_input.chr${chr}.bgen

        if test -f phewas_input.chr${chr}.bgen; then
            ${params.my_bgenix} -g phewas_input.chr${chr}.bgen -index -clobber
        else
            echo "${chr} returns no variants Please check snp list to match IDs on chromosome ${chr} or remove from chromosome list"
        fi

        ${params.my_bgenix} -g phewas_input.chr${chr}.bgen -index -clobber

        """
    stub:
        """
        touch phewas_input.chr${chr}.bgen
        touch phewas_input.chr${chr}.bgen.bgi
        """
}

process call_saige_step2_PLINK_binary {
    errorStrategy 'retry'
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --is_output_moreDetails=TRUE \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=0.05 \
         --LOCO=${params.LOCO} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
           > ${cohort_dir}.${pheno}.${chr}.txt

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt

        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bedFile=${plink_bed} \
         --bimFile=${plink_bim} \
         --famFile=${plink_fam} \
         --chrom=${chr} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_binary {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        """
        # Extract chromosome information
        ${params.my_bgenix} -g ${bgenFile} -i ${bgenFileIndex} -incl-range -list | head | tail -n +3 | awk -F"\t" 'NR==3{print \$3}' > chrs.txt

        file='chrs.txt'
        f1=0
        f2=0
        contains_chr=false
        contains_leading_zeros=false

        while IFS= read -r line; do
            if [[ "\$line" =~ chr ]]; then
                contains_chr=true
                f1=1
            fi
            if [[ "\$line" =~ 0[0-9] ]]; then
                contains_leading_zeros=true
                f2=1
            fi
        done < "\$file"

        if [[ \$f1 -ne 1 ]]; then
            contains_chr=false
        fi
        if [[ \$f2 -ne 1 ]]; then
            contains_leading_zeros=false
        fi

        if [[ \$contains_chr == true ]]; then
            chrprefix="chr"
        else
            chrprefix=''
        fi

        if [[ \$contains_leading_zeros == true && ${chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        # Call the R script with the evaluated chrom variable
        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bgenFile=${bgenFile} \
         --bgenFileIndex=${bgenFileIndex} \
         --AlleleOrder=ref-first \
         --sampleFile=${bgen_sample_file} \
         --chrom=\$CHROM_VAR \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
         > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file

    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        """
        # Extract chromosome information
        ${params.my_bgenix} -g ${bgenFile} -i ${bgenFileIndex} -incl-range -list | head | tail -n +3 | awk -F"\t" 'NR==3{print \$3}' > chrs.txt

        file='chrs.txt'
        f1=0
        f2=0
        contains_chr=false
        contains_leading_zeros=false

        while IFS= read -r line; do
            if [[ "\$line" =~ chr ]]; then
                contains_chr=true
                f1=1
            fi
            if [[ "\$line" =~ 0[0-9] ]]; then
                contains_leading_zeros=true
                f2=1
            fi
        done < "\$file"

        if [[ \$f1 -ne 1 ]]; then
            contains_chr=false
        fi
        if [[ \$f2 -ne 1 ]]; then
            contains_leading_zeros=false
        fi

        if [[ \$contains_chr == true ]]; then
            chrprefix="chr"
        else
            chrprefix=''
        fi

        if [[ \$contains_leading_zeros == true && ${chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --bgenFile=${bgenFile}    \
         --bgenFileIndex=${bgenFileIndex} \
         --sampleFile=${bgen_sample_file} \
         --AlleleOrder=ref-first \
         --chrom=\$CHROM_VAR \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt

        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_binary_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        //inputs
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
        """
              # Extract chromosome information
        ${params.my_bgenix} -g ${bgenFile} -i ${bgenFileIndex} -incl-range -list | head | tail -n +3 | awk -F"\t" 'NR==3{print \$3}' > chrs.txt

        file='chrs.txt'
        f1=0
        f2=0
        contains_chr=false
        contains_leading_zeros=false

        while IFS= read -r line; do
            if [[ "\$line" =~ chr ]]; then
                contains_chr=true
                f1=1
            fi
            if [[ "\$line" =~ 0[0-9] ]]; then
                contains_leading_zeros=true
                f2=1
            fi
        done < "\$file"

        if [[ \$f1 -ne 1 ]]; then
            contains_chr=false
        fi
        if [[ \$f2 -ne 1 ]]; then
            contains_leading_zeros=false
        fi

        if [[ \$contains_chr == true ]]; then
            chrprefix="chr"
        else
            chrprefix=''
        fi

        if [[ \$contains_leading_zeros == true && ${chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bgenFile=${bgenFile} \
         --bgenFileIndex=${bgenFileIndex} \
         --AlleleOrder=ref-first \
         --sampleFile=${bgen_sample_file}
         --chrom=\$CHROM_VAR \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=${params.firth_cutoff} \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """

    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}_${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)

        // actual inputs
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file
        tuple path(sparse_grm), path(sparse_grm_samples)

    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    shell:
    """
        # Extract chromosome information
        ${params.my_bgenix} -g ${bgenFile} -i ${bgenFileIndex} -incl-range -list | head | tail -n +3 | awk -F"\t" 'NR==3{print \$3}' > chrs.txt

        file='chrs.txt'
        f1=0
        f2=0
        contains_chr=false
        contains_leading_zeros=false

        while IFS= read -r line; do
            if [[ "\$line" =~ chr ]]; then
                contains_chr=true
                f1=1
            fi
            if [[ "\$line" =~ 0[0-9] ]]; then
                contains_leading_zeros=true
                f2=1
            fi
        done < "\$file"

        if [[ \$f1 -ne 1 ]]; then
            contains_chr=false
        fi
        if [[ \$f2 -ne 1 ]]; then
            contains_leading_zeros=false
        fi

        if [[ \$contains_chr == true ]]; then
            chrprefix="chr"
        else
            chrprefix=''
        fi

        if [[ \$contains_leading_zeros == true && ${chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 Rscript ${params.step2_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --bgenFile=${bgenFile}    \
         --bgenFileIndex=${bgenFileIndex} \
         --chrom=\$CHROM_VAR \
         --minMAF=${params.min_maf} \
         --minMAC=${params.min_mac} \
         --GMMATmodelFile=${step1_rda} \
         --varianceRatioFile=${step1_var} \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
           > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """

    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_binary_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
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
         --is_Firth_beta=TRUE \
         --pCutoffforFirth=0.05 \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
           > ${cohort_dir}.${pheno}.${chr}.log

         gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
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
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
         --LOCO=${params.LOCO} \
         --is_output_moreDetails=TRUE \
          > ${cohort_dir}.${pheno}.${chr}.log

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt.gz
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

include {
    paramToList
} from '../processes/saige_helpers.nf'

workflow SAIGE_VAR_STEP2 {
    take:
        step1_bin_output
        step1_quant_output
        use_step2_prefix
        bgen_sample_file
        IS_PHEWAS
    main:
        cohort = Channel.fromList(params.cohort_list)
        bin_pheno = Channel.fromList(paramToList(params.bin_pheno_list))
        quant_pheno = Channel.fromList(paramToList(params.quant_pheno_list))
        chromosome = Channel.fromList(params.chromosome_list)

        ftype = params.ftype

        // Build the chromosome-separated input file paths
        if (ftype == 'PLINK') {
            suffixes_list = ['.bed', '.bim', '.fam']
        } else if (ftype == 'BGEN') {
            suffixes_list = ['.bgen', '.bgen.bgi']
        }

        // File pattern = prefix + chr + ext
        step2_chr_sep_input = chromosome.map { chr -> \
            new Tuple(chr, new Tuple(*suffixes_list.collect {
                ext -> "${use_step2_prefix}${chr}${ext}" })) }
        //step2_chr_sep_input.view{"s2 chr: ${it}"}
        if (IS_PHEWAS) {
            // with PheWAS, we first filter the input files to significantly speed up SAIGE run time
            if (ftype == 'PLINK') {
                // replace with filtered plink files
                step2_chr_sep_input = filter_snps_plink(step2_chr_sep_input, params.snplist)
            } else if (ftype == 'BGEN') {
                // replace with filtered bgen files
                step2_chr_sep_input = filter_snps_bgen(step2_chr_sep_input, params.snplist)
            }
        }

        // For the join, we need to combine our group files with cohort and phenotype
        // These will represent all phenos x cohorts x chromosomes
        geno_data_parallel_bin = step2_chr_sep_input.combine(cohort).combine(bin_pheno).map {
            chr, geno_fileset, cohort, pheno -> new Tuple(cohort, pheno, chr, geno_fileset)
        }
        geno_data_parallel_quant = step2_chr_sep_input.combine(cohort).combine(quant_pheno).map {
            chr, geno_fileset, cohort, pheno -> new Tuple(cohort, pheno, chr, geno_fileset)
        }

        // These have KEPT (phenos, cohorts) x chromosomes
        step2_bin_input = step1_bin_output.combine(chromosome)
        step2_quant_input = step1_quant_output.combine(chromosome)

        // These sections synchronizes our genetic data input files to our step2 input by joining on (cohort, pheno, chr)
        // We're using the join as a filter essentially
        chr_geno_files_bin = step2_bin_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(geno_data_parallel_bin, by: [0, 1, 2]).map {
            cohort, pheno, chr, geno_fileset -> geno_fileset
        }

        // Same for quantitative phenotypes
        chr_geno_files_quant = step2_quant_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(geno_data_parallel_quant, by: [0, 1, 2]).map {
            cohort, pheno, chr, geno_fileset -> geno_fileset
        }

        if (params.use_sparse_GRM) {
            // Define sparse GRM input
            sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)
            if (ftype == 'PLINK') {
                // Call plink sparse GRM processes
                // Binary phenotypes + PLINK input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK_binary_with_sparse_GRM(
                    step2_bin_input, chr_geno_files_bin, sparse_grm_input)
                // Quantitative phenotypes + PLINK input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK_quant_with_sparse_GRM(
                    step2_quant_input, chr_geno_files_quant, sparse_grm_input)
            } else if (ftype == 'BGEN') {
                // Call bgen sparse GRM processes
                // Binary phenotypes + BGEN input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_BGEN_binary_with_sparse_GRM(
                    step2_bin_input, chr_geno_files_bin, bgen_sample_file, sparse_grm_input)
                // Quantitative phenotypes + BGEN input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_BGEN_quant_with_sparse_GRM(
                    step2_quant_input, chr_geno_files_quant, bgen_sample_file, sparse_grm_input)
            }
        } else {
            // Don't use GRM input
            if (ftype == 'PLINK') {
                // Call plink regular processes
                // Binary phenotypes + PLINK input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK_binary(
                    step2_bin_input, chr_geno_files_bin)
                // Quantitative phenotypes + PLINK input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK_quant(
                    step2_quant_input, chr_geno_files_quant)
            } else if (ftype == 'BGEN') {
                // Call bgen regular processes
                // Binary phenotypes + BGEN input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_BGEN_binary(
                    step2_bin_input, chr_geno_files_bin, bgen_sample_file)
                // Quantitative phenotypes + BGEN input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_BGEN_quant(
                    step2_quant_input, chr_geno_files_quant, bgen_sample_file)
            }
        }
    emit:
        step2_bin_output
        step2_quant_output
        }
