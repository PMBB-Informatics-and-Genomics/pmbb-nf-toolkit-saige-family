/*
This is a collection of SAIGE Step 2 processes
*/

include {
    paramToList
    readChunksManifest
    WriteSaigeStep2Command
} from '../processes/saige_helpers.nf'

process filter_snps_plink {
    publishDir "${launchDir}/Filtered_Files/"
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
            { echo "${chr} returns no variants Please check snp list to match IDs on chromosome ${chr} or remove from chromosome list" ; exit ; }
            
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
    label 'saige_process'
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --bedFile=${plink_bed} \
            --bimFile=${plink_bim} \
            --famFile=${plink_fam} \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

            gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
            """
        stub:
            """
            touch ${cohort_dir}.${pheno}.${chr}.txt
            gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
            touch ${cohort_dir}.${pheno}.${chr}.log
            """
}

process call_saige_step2_PLINK_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --chrom=\$CHROM_VAR \
        --GMMATmodelFile=${step1_rda} \
        --varianceRatioFile=${step1_var} \
        --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        --minMAF=${params.min_maf} \
        --minMAC=${params.min_mac} \
        --LOCO=${params.LOCO} \
        --is_imputed_data=${params.isImputed} \
        --minInfo=${params.minInfo} \
        --is_output_moreDetails=TRUE \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_survival {
    errorStrategy 'retry'
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash

        set -euo pipefail

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --bedFile=${plink_bed} \
            --bimFile=${plink_bim} \
            --famFile=${plink_fam} \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

            gzip -9 ${cohort_dir}.${pheno}.${chr}.txt

            """
        stub:
            """
            touch ${cohort_dir}.${pheno}.${chr}.txt
            gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
            touch ${cohort_dir}.${pheno}.${chr}.log
            """
}

process call_saige_step2_BGEN_binary {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    memory = '100 GB'
    label 'saige_process'
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }

    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        # Call the R script with the evaluated chrom variable
        stdbuf -e0 -o0 ${params.step2_script} \
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
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --is_output_moreDetails=TRUE \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    label 'saige_process'
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
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr

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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 ${params.step2_script} \
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
         --is_imputed_data=${params.isImputed} \
         --minInfo=${params.minInfo} \
         --is_output_moreDetails=TRUE \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt

        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_survival {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(bgenFile), path(bgenFileIndex)
        path bgen_sample_file
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        # Call the R script with the evaluated chrom variable
        stdbuf -e0 -o0 ${params.step2_script} \
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
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --is_output_moreDetails=TRUE \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_binary_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
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
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 ${params.step2_script} \
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
         --is_imputed_data=${params.isImputed} \
         --minInfo=${params.minInfo} \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """

    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}_${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
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
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 ${params.step2_script} \
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
         --is_imputed_data=${params.isImputed} \
         --minInfo=${params.minInfo} \
         --is_output_moreDetails=TRUE \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """

    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_BGEN_survival_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = bgenFile.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
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
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr

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

        if [[ \$contains_leading_zeros == true && ${stripped_chr} -lt 10 ]]; then
            iszero=0
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt

        export CHROM_VAR=\$chrom

        stdbuf -e0 -o0 ${params.step2_script} \
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
         --is_imputed_data=${params.isImputed} \
         --minInfo=${params.minInfo} \
         --is_output_moreDetails=TRUE \
         --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """

    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}_${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_binary_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr

        """
        #!/bin/bash

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
        --sparseGRMFile=${sparse_grm} \
        --sparseGRMSampleIDFile=${sparse_grm_samples} \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --chrom=\$CHROM_VAR \
        --minMAF=${params.min_maf} \
        --minMAC=${params.min_mac} \
        --GMMATmodelFile=${step1_rda} \
        --varianceRatioFile=${step1_var} \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=${params.LOCO} \
        --is_imputed_data=${params.isImputed} \
        --minInfo=${params.minInfo} \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr

        """
        #!/bin/bash

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
        --sparseGRMFile=${sparse_grm} \
        --sparseGRMSampleIDFile=${sparse_grm_samples} \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --chrom=\$CHROM_VAR \
        --minMAF=${params.min_maf} \
        --minMAC=${params.min_mac} \
        --GMMATmodelFile=${step1_rda} \
        --varianceRatioFile=${step1_var} \
        --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        --LOCO=${params.LOCO} \
        --is_imputed_data=${params.isImputed} \
        --minInfo=${params.minInfo} \
        --is_output_moreDetails=TRUE \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK_survival_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = plink_bed.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        // variables
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(plink_bed), path(plink_bim), path(plink_fam)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash

        # Input file passed as an argument
        
        # Read the first line and extract the first column
        first_column=\$(head -n1 ${plink_bim} | cut -f1)

        # Initialize variables
        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""
        

        # Check for "chr"
        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        # Check for leading zeros and numeric values
        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        # Determine the chromosome prefix and leading zero
        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        # Extract the numeric part of the chromosome, assuming it follows "chr" (if applicable)
        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"  # Extract only the digits
        else
            chr="0"  # Default value if no valid chromosome number is found
        fi

        # Check if the chromosome is less than 10 and leading zeros are present
        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        # Construct the final chromosome string
        chrom="\${chrprefix}\${iszero}${stripped_chr}"

        # Output the chromosome string to a file
        echo "\$chrom" > chrom.txt

        # Export the chromosome variable for Nextflow
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
        --sparseGRMFile=${sparse_grm} \
        --sparseGRMSampleIDFile=${sparse_grm_samples} \
        --bedFile=${plink_bed} \
        --bimFile=${plink_bim} \
        --famFile=${plink_fam} \
        --chrom=\$CHROM_VAR \
        --minMAF=${params.min_maf} \
        --minMAC=${params.min_mac} \
        --GMMATmodelFile=${step1_rda} \
        --varianceRatioFile=${step1_var} \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=${params.LOCO} \
        --is_imputed_data=${params.isImputed} \
        --minInfo=${params.minInfo} \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
        > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}


process call_saige_step2_PLINK2_binary {
    errorStrategy 'retry'
    label 'saige_process'
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK2_quant {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --is_output_moreDetails=TRUE \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK2_survival {
    errorStrategy 'retry'
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK2_binary_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --sparseGRMFile=${sparse_grm} \
            --sparseGRMSampleIDFile=${sparse_grm_samples} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK2_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --sparseGRMFile=${sparse_grm} \
            --sparseGRMSampleIDFile=${sparse_grm_samples} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --is_output_moreDetails=TRUE \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}

process call_saige_step2_PLINK2_survival_with_sparse_GRM {
    publishDir "${launchDir}/${cohort_dir}/Saige_Step2_Results/"
    label 'saige_process'
    disk {
        def fileSizeGb = pgen_file.size() / (1024 ** 3)
        return "${50 + fileSizeGb} GB"
    }
    input:
        tuple val(cohort_dir), val(pheno), path(step1_rda), path(step1_var), val(chr)
        tuple path(pgen_file), path(pvar_file), path(psam_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.txt.gz")
        tuple val(cohort_dir), val(pheno), val(chr), path("${cohort_dir}.${pheno}.${chr}.log")
    script:
        def stripped_chr = chr.contains('_') ? chr.toString().split('_')[0] : chr
        """
        #!/bin/bash
        set -euo pipefail

        first_column=\$(grep -vm1 '^#' ${pvar_file} | cut -f1)

        contains_chr=false
        contains_leading_zeros=false
        chrprefix=""
        iszero=""

        if [[ "\$first_column" == *"chr"* ]]; then
            contains_chr=true
        fi

        if [[ "\$first_column" =~ ^0[0-9]+ ]]; then
            contains_leading_zeros=true
        fi

        if [[ "\$contains_chr" == true ]]; then
            chrprefix="chr"
        else
            chrprefix=""
        fi

        if [[ "\$first_column" =~ [0-9]+ ]]; then
            chr="\${first_column//[^0-9]/}"
        else
            chr="0"
        fi

        if [[ "\$contains_leading_zeros" == true && ${stripped_chr} -lt 10 ]]; then
            iszero="0"
        else
            iszero=""
        fi

        chrom="\${chrprefix}\${iszero}${stripped_chr}"
        echo "\$chrom" > chrom.txt
        export CHROM_VAR="\$chrom"

        stdbuf -e0 -o0 ${params.step2_script} \
            --sparseGRMFile=${sparse_grm} \
            --sparseGRMSampleIDFile=${sparse_grm_samples} \
            --pgenPrefix=${pgen_file.baseName} \
            --AlleleOrder=ref-first \
            --chrom=\$CHROM_VAR \
            --is_output_moreDetails=TRUE \
            --minMAF=${params.min_maf} \
            --minMAC=${params.min_mac} \
            --GMMATmodelFile=${step1_rda} \
            --varianceRatioFile=${step1_var} \
            --is_Firth_beta=TRUE \
            --pCutoffforFirth=0.05 \
            --LOCO=${params.LOCO} \
            --is_imputed_data=${params.isImputed} \
            --minInfo=${params.minInfo} \
            --SAIGEOutputFile=${cohort_dir}.${pheno}.${chr}.txt \
            > ${cohort_dir}.${pheno}.${chr}.log 2>&1

        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        """
    stub:
        """
        touch ${cohort_dir}.${pheno}.${chr}.txt
        gzip -9 ${cohort_dir}.${pheno}.${chr}.txt
        touch ${cohort_dir}.${pheno}.${chr}.log
        """
}


workflow SAIGE_VAR_STEP2 {
    take:
        step1_bin_output
        step1_quant_output
        step1_survival_output
        step2_genetic_data_prefix
        bgen_sample_file
        IS_PHEWAS
    main:
        params.isImputed = "FALSE"  // Default to FALSE if not provided
        params.minInfo = 0      // Default to 0.3 if not provided
        cohort = Channel.fromList(params.cohort_list)
        bin_pheno = Channel.fromList(paramToList(params.bin_pheno_list))
        quant_pheno = Channel.fromList(paramToList(params.quant_pheno_list))
        survival_pheno = Channel.fromList(paramToList(params.survival_pheno_list))
        ftype = params.ftype

        // Build the chromosome-separated input file paths
        if (ftype == 'PLINK1') {
            suffixes_list = ['.bed', '.bim', '.fam']
        } else if (ftype == 'PLINK2') {
            suffixes_list = ['.pgen', '.pvar', '.psam']
        } else if (ftype == 'BGEN') {
            suffixes_list = ['.bgen', '.bgen.bgi']
        }

        println("Step 2 has been triggered. Input channels being created")

        if (params.enable_chunking) {
            // Manifest provides [chrom, file_stem] rows; chromosome column is already stripped of 'chr' prefix
            def all_manifest = readChunksManifest(params.chunks_manifest.toString())
            chromosome_set = params.chromosome_list.collect { it.toString() }.toSet()
            filtered_manifest = all_manifest.findAll { row -> chromosome_set.contains(row[0]) }
            chromosome = Channel.fromList(filtered_manifest)  // emits [chrom, stem] lists

            // Build file tuples directly from the manifest stem; use chrom + basename as unique chunk key
            step2_chr_sep_input = chromosome.map { row ->
                def chrom = row[0]
                def stem = row[1]
                def chunk_key = "${chrom}_${stem.tokenize('/').last()}"
                new Tuple(chunk_key, new Tuple(*suffixes_list.collect { ext -> "${stem}${ext}" }))
            }
        } else {
            chromosome = Channel.fromList(params.chromosome_list)

            // File pattern = prefix + chr + ext
            step2_chr_sep_input = chromosome.map { chr -> \
                new Tuple(chr, new Tuple(*suffixes_list.collect {
                    ext -> "${step2_genetic_data_prefix}${chr}${ext}" })) }
        }
        // expected results: [chr, [file1, file2, ...]]
        if (IS_PHEWAS) {
            // with PheWAS, we first filter the input files to significantly speed up SAIGE run time
            if (ftype == 'PLINK1') {
                // replace with filtered plink files
                step2_chr_sep_input = filter_snps_plink(step2_chr_sep_input, params.snplist)
            } else if (ftype == 'PLINK2') {
                throw new Exception('filter_snps for PLINK2 is not yet supported in PheWAS workflows')
            } else if (ftype == 'BGEN') {
                // replace with filtered bgen files
                step2_chr_sep_input = filter_snps_bgen(step2_chr_sep_input, params.snplist)
            }
        }

        // For the join, we need to combine our group files with cohort and phenotype
        // These will represent all phenos x cohorts x chromosomes
        // Expected: [chr, [file1, file2, ...]] -> [cohort, pheno, chr, geno_fileset]
        geno_data_parallel_bin = step2_chr_sep_input.combine(cohort).combine(bin_pheno).map {
            chr, geno_fileset, cohort, pheno -> new Tuple(cohort, pheno, chr, geno_fileset)
        }
        geno_data_parallel_quant = step2_chr_sep_input.combine(cohort).combine(quant_pheno).map {
            chr, geno_fileset, cohort, pheno -> new Tuple(cohort, pheno, chr, geno_fileset)
        }
        geno_data_parallel_survival = step2_chr_sep_input.combine(cohort).combine(survival_pheno).map {
            chr, geno_fileset, cohort, pheno -> new Tuple(cohort, pheno, chr, geno_fileset)
        }

        // Expected: [cohort, pheno, rda, var] -> [cohort, pheno, rda, var, chr]
        // Use chr_keys (not chromosome) so chunking mode doesn't inject the stem as an extra element
        chr_keys = step2_chr_sep_input.map { chr, geno_fileset -> chr }
        step2_bin_input = step1_bin_output.combine(chr_keys)
        step2_quant_input = step1_quant_output.combine(chr_keys)
        step2_survival_input = step1_survival_output.combine(chr_keys)

        println("Input channels created")
        
        // These sections synchronizes our genetic data input files to our step2 input by joining on (cohort, pheno, chr)
        // expected: [cohort, pheno, rda, var, chr] + [cohort, pheno, chr, geno_fileset] -> [cohort, pheno, chr, geno_fileset] -> geno_fileset

        // For binary phenotypes
        chr_geno_files_bin = step2_bin_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(geno_data_parallel_bin, by: [0, 1, 2]).map {
            cohort, pheno, chr, geno_fileset -> geno_fileset
        }

        // For quantitative phenotypes
        chr_geno_files_quant = step2_quant_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(geno_data_parallel_quant, by: [0, 1, 2]).map {
            cohort, pheno, chr, geno_fileset -> geno_fileset
        }

        // For survival phenotypes
        chr_geno_files_survival = step2_survival_input.map {
            cohort, pheno, rda, var, chr -> new Tuple(cohort, pheno, chr)
        }.join(geno_data_parallel_survival, by: [0, 1, 2]).map {
            cohort, pheno, chr, geno_fileset -> geno_fileset
        }

        println("Step 2 Process being called.")

        if (params.use_sparse_GRM) {
            // Define sparse GRM input
            sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)
            if (ftype == 'PLINK1') {
                // Call plink sparse GRM processes
                // Binary phenotypes + PLINK1 input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK_binary_with_sparse_GRM(
                    step2_bin_input, chr_geno_files_bin, sparse_grm_input)
                // Quantitative phenotypes + PLINK1 input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK_quant_with_sparse_GRM(
                    step2_quant_input, chr_geno_files_quant, sparse_grm_input)
                // Survival phenotypes + PLINK1 input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_PLINK_survival_with_sparse_GRM(
                    step2_survival_input, chr_geno_files_survival, sparse_grm_input)

            } else if (ftype == 'PLINK2') {
                // Call plink2 sparse GRM processes
                // Binary phenotypes + PLINK2 input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK2_binary_with_sparse_GRM(
                    step2_bin_input, chr_geno_files_bin, sparse_grm_input)
                // Quantitative phenotypes + PLINK2 input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK2_quant_with_sparse_GRM(
                    step2_quant_input, chr_geno_files_quant, sparse_grm_input)
                // Survival phenotypes + PLINK2 input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_PLINK2_survival_with_sparse_GRM(
                    step2_survival_input, chr_geno_files_survival, sparse_grm_input)

            } else if (ftype == 'BGEN') {
                // Call bgen sparse GRM processes
                // Binary phenotypes + BGEN input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_BGEN_binary_with_sparse_GRM(
                    step2_bin_input, chr_geno_files_bin, bgen_sample_file, sparse_grm_input)
                // Quantitative phenotypes + BGEN input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_BGEN_quant_with_sparse_GRM(
                    step2_quant_input, chr_geno_files_quant, bgen_sample_file, sparse_grm_input)
                // Survival phenotypes + BGEN input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_BGEN_survival_with_sparse_GRM(
                    step2_survival_input, chr_geno_files_survival, bgen_sample_file, sparse_grm_input)
                
                
            }
        } else {
            // Don't use GRM input
            if (ftype == 'PLINK1') {
                // Call plink regular processes
                // Binary phenotypes + PLINK1 input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK_binary(
                    step2_bin_input, chr_geno_files_bin)
                // Quantitative phenotypes + PLINK1 input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK_quant(
                    step2_quant_input, chr_geno_files_quant)
                // Survival phenotypes + PLINK1 input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_PLINK_survival(
                    step2_survival_input, chr_geno_files_survival)
            } else if (ftype == 'PLINK2') {
                // Call plink2 regular processes
                // Binary phenotypes + PLINK2 input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_PLINK2_binary(
                    step2_bin_input, chr_geno_files_bin)
                // Quantitative phenotypes + PLINK2 input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_PLINK2_quant(
                    step2_quant_input, chr_geno_files_quant)
                // Survival phenotypes + PLINK2 input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_PLINK2_survival(
                    step2_survival_input, chr_geno_files_survival)
            } else if (ftype == 'BGEN') {
                // Call bgen regular processes
                // Binary phenotypes + BGEN input
                (step2_bin_output, step2_bin_logs) = call_saige_step2_BGEN_binary(
                    step2_bin_input, chr_geno_files_bin, bgen_sample_file)
                // Quantitative phenotypes + BGEN input
                (step2_quant_output, step2_bin_logs) = call_saige_step2_BGEN_quant(
                    step2_quant_input, chr_geno_files_quant, bgen_sample_file)
                // Survival phenotypes + BGEN input
                (step2_survival_output, step2_survival_logs) = call_saige_step2_BGEN_survival(
                    step2_survival_input, chr_geno_files_survival, bgen_sample_file)
            }
        }
    emit:
        step2_bin_output
        step2_quant_output
        step2_survival_output
        }
