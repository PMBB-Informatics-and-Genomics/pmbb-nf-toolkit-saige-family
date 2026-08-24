List<List<String>> readChunksManifest(String manifestPath) {
    def f = new File(manifestPath)
    if (!f.exists()) {
        println("ERROR: chunks_manifest file not found: ${manifestPath}")
        System.exit(1)
    }
    def parentDir = f.parentFile.absolutePath
    return f.readLines()
        .findAll { line -> line.trim() && !line.startsWith('#') }
        .collect { line ->
            def parts = line.split('\t')
            if (parts.size() < 2) {
                println("ERROR: Invalid chunks_manifest line (expected chromosome<TAB>file_path_stem): ${line}")
                System.exit(1)
            }
            def chrom = parts[0].trim().replace('chr', '')
            def stem = parts[1].trim()
            // Resolve paths relative to the manifest file's directory
            if (!stem.startsWith('/')) {
                stem = "${parentDir}/${stem}"
            }
            return [chrom, stem]
        }
}

List paramToList(param) {
    // This function takes a parameter and returns a list
    // If the parameter is already a list, it returns it as is
    if (param instanceof List) {
        return param
    // Handle null cases
    } else if (param == '') {
        return []
    } else if (param == null) {
        return []
    } else {
        // If the parameter is a string, check if it is a file
        def file = new File(param.toString())
        if (file.exists()) {
            return file.readLines()
        } else {
            return [param.toString()]
        }
    }
}

Map get_script_file_names() {
    // This method kind of serves as a manifest for all of our scripts rather than completely hard-coding the paths

    script_names = [:]
    script_names['cohort_setup'] = "${moduleDir}/../scripts/set_up_cohort_directory.py"
    script_names['pheno_table'] = "${moduleDir}/../scripts/make_pheno_summary_table.py"
    script_names['pheno_covar_plots'] = "${moduleDir}/../scripts/make_pheno_covar_summary_plots.py"
    script_names['exwas_methods'] = "${moduleDir}/../scripts/generate_exwas_methods.py"
    script_names['exwas_singles_plots'] = "${moduleDir}/../scripts/make_saige_exwas_singles_plots.py"
    script_names['exwas_regions_plots'] = "${moduleDir}/../scripts/make_saige_exwas_regions_plots.py"
    script_names['exwas_report'] = "${moduleDir}/../scripts/make_exwas_report.py"
    script_names['gwas_report'] = "${moduleDir}/../scripts/generate_gwas_html_report.py"
    script_names['gwas_plots_with_annot'] = "${moduleDir}/../scripts/make_saige_gwas_plots_annotate.py"
    script_names['gwas_plots'] = "${moduleDir}/../scripts/make_saige_gwas_singles_plots.py"
    script_names['merge'] = "${moduleDir}/../scripts/merge_and_filter_saige_results.py"
    script_names['gb_phewas_plots'] = "${moduleDir}/../scripts/make_saige_gene_burden_phewas_plots_v2.py"
    script_names['v_phewas_plots'] = "${moduleDir}/../scripts/make_saige_var_phewas_plots.py"
    return script_names
}

Map check_input_genetic_data_parameters(params, pipeline) {

    genetic_data_params_clean = [:]

    if (params.use_sparse_GRM) {
        // Required args when using the sparse GRM are:
        // step1_plink_prefix, step1_sparse_grm, step1_sparse_grm_samples
        if (params.step1_plink_prefix == null) {
            println('You set use_sparse_GRM = true, which means that step1_plink_prefix MUST be set.')
            System.exit(1)
        }
        if (params.step1_sparse_grm == null) {
            println('You set use_sparse_GRM = true, which means that step1_sparse_grm MUST be set.')
            System.exit(1)
        }
        if (params.step1_sparse_grm_samples == null) {
            println('You set use_sparse_GRM = true, which means that step1_sparse_grm_samples MUST be set.')
            System.exit(1)
        }
        genetic_data_params_clean.use_step1_prefix = params.step1_plink_prefix
    } else {
        if (pipeline == 'GWAS' || pipeline == 'Variant PheWAS') {
            // Required argument for step 1 without the GRM is step1_plink_prefix
            if (params.step1_plink_prefix == null) {
                println("You are running a ${pipeline} without the sparse GRM, which means that step1_plink_prefix MUST be set.")
                System.exit(1)
            }
            genetic_data_params_clean.use_step1_prefix = params.step1_plink_prefix
        } else {
            // For ExWAS and Gene PheWAS things are a little more tricky
            // Given that no sparse GRM is used,
            // people can provide JUST exome OR (step1_plink_prefix AND step2_plink_prefix)
            if (params.exome_plink_prefix != null && (params.step1_plink_prefix == null && params.step2_plink_prefix == null)) {
                genetic_data_params_clean.use_step1_prefix = params.exome_plink_prefix
            } else if (params.exome_plink_prefix == null && (params.step1_plink_prefix != null && params.step2_plink_prefix != null)) {
                genetic_data_params_clean.use_step1_prefix = params.step1_plink_prefix
            } else {
                println("You are running a(n) ${pipeline} without the sparse GRM, which means that: \
                 either exome_plink_prefix OR (step1_plink_prefix AND step2_plink_prefix) MUST be set.")
                System.exit(1)
            }
        }
    }

    // Now deal with step 2
    if (pipeline == 'GWAS' || pipeline == 'Variant PheWAS') {
        // Required arguments for step 2 are:
        // step2_plink_prefix, step2_bgen_prefix, bgen_samplefile
        // depending on ftype
        genetic_data_params_clean.is_chr_separated = true

        if (params.ftype == 'PLINK1') {
            if (params.step2_plink_prefix == null) {
                println("You are running a(n) ${pipeline} with PLINK1 file type, which means that step2_plink_prefix (chr-separated) MUST be set.")
                System.exit(1)
            }
            genetic_data_params_clean.use_step2_prefix = params.step2_plink_prefix
        } else if (params.ftype == 'PLINK2') {
            if (params.step2_pgen_prefix == null) {
                println("You are running a(n) ${pipeline} with PLINK2 file type, which means that step2_pgen_prefix (chr-separated) MUST be set.")
                System.exit(1)
            }
            genetic_data_params_clean.use_step2_prefix = params.step2_pgen_prefix
        } else if (params.ftype == 'BGEN') {
            if (params.step2_bgen_prefix == null) {
                println("You are running a(n) ${pipeline} with BGEN file type, which means that step2_bgen_prefix (chr-separated) MUST be set.")
                System.exit(1)
            }
            if (params.bgen_samplefile == null) {
                println("You are running a(n) ${pipeline} with BGEN file type, which means that bgen_samplefile MUST be set.")
                System.exit(1)
            }
            genetic_data_params_clean.use_step2_prefix = params.step2_bgen_prefix
        } else {
            println("params.ftype must be either 'PLINK1', 'PLINK2', or 'BGEN'")
            System.exit(1)
        }
    } else {
        // For ExWAS and Gene PheWAS things are a little more tricky
        // People can provide JUST exome OR (step1_plink_prefix AND step2_plink_prefix)
        if (!params.use_sparse_GRM && params.exome_plink_prefix != null && (params.step1_plink_prefix == null && params.step2_plink_prefix == null)) {
            // No Sparse GRM, exome prefix only
            genetic_data_params_clean.use_step2_prefix = params.exome_plink_prefix
            genetic_data_params_clean.is_chr_separated = false
        } else if (!params.use_sparse_GRM && params.exome_plink_prefix == null && (params.step2_plink_prefix != null)) {
            // No Sparse GRM, step1 and step2 prefix
            genetic_data_params_clean.use_step2_prefix = params.step2_plink_prefix
            genetic_data_params_clean.is_chr_separated = true
        } else if (params.use_sparse_GRM && params.exome_plink_prefix != null && params.step2_plink_prefix == null) {
            // GRM, exome prefix not step2 prefix
            genetic_data_params_clean.use_step2_prefix = params.exome_plink_prefix
            genetic_data_params_clean.is_chr_separated = false
        } else if (params.use_sparse_GRM && params.exome_plink_prefix == null && params.step2_plink_prefix != null) {
            // GRM, step2 prefix not exome prefix
            genetic_data_params_clean.use_step2_prefix = params.step2_plink_prefix
            genetic_data_params_clean.is_chr_separated = true
        } else {
            println("You are running a(n) ${pipeline} which means that: either exome_plink_prefix OR (step1_plink_prefix AND step2_plink_prefix) MUST be set.")
            System.exit(1)
        }
    }

    return genetic_data_params_clean
}

def WriteSaigeStep2Command(ftype, file_set_prefix, bgen_sample_file, is_chr_separated, out_prefix, step2_script, step1_rda, step1_var, cohort_dir, pheno) {
    if (!(ftype in ['PLINK1','PLINK2','BGEN'])) {
        throw new IllegalArgumentException("ftype must be either 'PLINK1', 'PLINK2', or 'BGEN'")
    }

    // Build file-specific args using prefixes/placeholders
    def fileArgs = ''
    if (ftype == 'PLINK1') {
        if (is_chr_separated) {
            fileArgs = "--plinkFile ${file_set_prefix} --chromosomeNum \$CHROM_VAR"
        } else {
            fileArgs = "--plinkFile ${file_set_prefix}"
        }
    } else if (ftype == 'PLINK2') {
        if (is_chr_separated) {
            fileArgs = "--pgenFile ${file_set_prefix}.pgen --pvarFile ${file_set_prefix}.pvar --psamFile ${file_set_prefix}.psam --chrom \$CHROM_VAR"
        } else {
            fileArgs = "--pgenFile ${file_set_prefix}.pgen --pvarFile ${file_set_prefix}.pvar --psamFile ${file_set_prefix}.psam"
        }
    } else { // BGEN
        def bgenFile = file_set_prefix
        def bgenIndex = "${file_set_prefix}.bgi"
        if (is_chr_separated) {
            fileArgs = "--bgenFile=${bgenFile} --bgenFileIndex=${bgenIndex} --sampleFile=${bgen_sample_file} --chrom=\$CHROM_VAR"
        } else {
            fileArgs = "--bgenFile=${bgenFile} --bgenFileIndex=${bgenIndex} --sampleFile=${bgen_sample_file}"
        }
    }
    
    // Use literal $CHROM_VAR in output paths (will be evaluated by bash in process)
    def outputFile = "${cohort_dir}.${pheno}.\$CHROM_VAR.txt"
    def logFile = "${cohort_dir}.${pheno}.\$CHROM_VAR.log"
    
    def commonArgs = "--GMMATmodelFile=${step1_rda} --varianceRatioFile=${step1_var} --LOCO=\${params.LOCO} --is_output_moreDetails=TRUE ${additional_args}"
    
    def cmd = """\
    stdbuf -e0 -o0 Rscript ${step2_script} \\
        ${fileArgs} \\
        ${commonArgs} \\
        --SAIGEOutputFile=${outputFile} \\
        > ${logFile}
    gzip -9 ${outputFile}
    """.stripIndent()
    
    return cmd
}

def getGpuMachineTypeChannel(numVariants, numParticipants) {
    requiredMemoryGB = (4 * numVariants * numParticipants) / 1e9
    Map gpuMemoryToMachineType = [
    40: 'a2-highgpu-1g',    // 1 GPU, 40GB total GPU memory
    80: 'a2-highgpu-2g',    // 2 GPUs, 80GB total GPU memory
    160: 'a2-highgpu-4g',   // 4 GPUs, 160GB total GPU memory
    320: 'a2-highgpu-8g',   // 8 GPUs, 320GB total GPU memory
    640: 'a2-megagpu-16g'   // 16 GPUs, 640GB total GPU memory
    ]
    // Sorted list of memory options
    def memoryOptions = gpuMemoryToMachineType.keySet().sort()
    // Find the smallest option that meets or exceeds the requirement
    def selectedMemory = null
    for (memory in memoryOptions) {
        if (memory >= requiredMemoryGB) {
            selectedMemory = memory
            break
        }
    }   
    // If no option is large enough, print warning and exit
    if (selectedMemory == null) {
        def maxAvailable = memoryOptions.isEmpty() ? 0 : memoryOptions[-1]
        log.warn "ERROR: Required GPU memory (${requiredMemoryGB}GB) exceeds maximum available (${maxAvailable}GB)"
        error "Workflow aborted due to insufficient GPU memory"
    }
    // Grab machine type and return as a channel
    def machineType = gpuMemoryToMachineType[selectedMemory]
    return machineType
}

import groovy.json.JsonBuilder
process dump_params_to_json {
    publishDir "${launchDir}/Summary", mode: 'copy'

    input:
        val params_dict
        val pipeline_name
    output:
        path("${pipeline_name}_params.json")
    shell:
        """
        echo '${new JsonBuilder(params_dict).toPrettyString().replace(';', '|')}' > ${pipeline_name}_params.json
        """
}

