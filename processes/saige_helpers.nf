List paramToList(param) {
    // Function to check if parameter is a file or list
    if (param instanceof List) {
        // If the parameter is a list, return list
        return param
    } else if (param instanceof String && new File(param).exists()) {
        // If the parameter is a file, read the file contents into a list
        listFromFile = new File(param).readLines()
        return listFromFile
    } else if (param == null) {
        // If the parameter is null, return an empty list
        return []
    }
    // Handle cases where the parameter is neither a list nor a valid file
    // throw new IllegalArgumentException("Parameter must be a list or an existing file")
    return param
}

Map get_script_file_names() {
    // This method kind of serves as a manifest for all of our scripts rather than completely hard-coding the paths

    script_names = [:]
    script_names['cohort_setup'] = "${projectDir}/../scripts/set_up_cohort_directory.py"
    script_names['pheno_table'] = "${projectDir}/../scripts/make_pheno_summary_table.py"
    script_names['pheno_covar_plots'] = "${projectDir}/../scripts/make_pheno_covar_summary_plots.py"
    script_names['exwas_methods'] = "${projectDir}/../scripts/generate_exwas_methods.py"
    script_names['exwas_singles_plots'] = "${projectDir}/../scripts/make_saige_exwas_singles_plots.py"
    script_names['exwas_regions_plots'] = "${projectDir}/../scripts/make_saige_exwas_regions_plots.py"
    script_names['gwas_plots_with_annot'] = "${projectDir}/../scripts/make_saige_gwas_plots_annotate.py"
    script_names['gwas_plots'] = "${projectDir}/../scripts/make_saige_gwas_plots.py"
    script_names['merge'] = "${projectDir}/../scripts/merge_and_filter_saige_results.py"

    return script_names
}
