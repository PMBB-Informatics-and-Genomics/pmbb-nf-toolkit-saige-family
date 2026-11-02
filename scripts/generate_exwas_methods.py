import json
import dominate
from dominate.tags import *

def addCssJsToHead():
    '''
    Adds bootstrap CSS resources to the HTML header so we can use their styling elements
    '''
    link(href='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css',
        rel='stylesheet',
        integrity='sha384-4bw+/aepP/YC94hEpVNVgiZdgIC5+VKNBQNGCHeKRQN+PtmoHDEXuppvnDJzQIu9',
        crossorigin='anonymous')
    script(src='https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/js/bootstrap.bundle.min.js',
            integrity='sha384-HwwvtgBNo3bZJJLYd8oVXjrBZt8cqVSpeBNS5n7C8IVInixGAoxmnlMuBnhbgrkm',
            crossorigin='anonymous')
    meta(charset='utf-8', name='viewport', content='width=device-width')

# Get dictionary from JSON file generated in the process pre-script
params = json.load(open('pipeline_params.txt'))
print(params)

# Start HTML document
outfile = dominate.document(title='ExWAS Methods Summary')

with outfile.head:
    # Header
    addCssJsToHead()

# Start nested HTML Elements
with outfile:
    h1('ExWAS Methods Summary')

    # Software versions section
    with div():
        h2('Software Versions')
        p('SAIGE Version: wzhou88/saige:1.2.0')
        p(f'Python .exe: {params["my_python"]}')
        p(f'Plink: plink/2.0-20210505')
    
    # Workflow-wide parameters
    with div():
        h2('Workflow Set-Up')
        p(f'Cohorts: {",".join(params["cohort_list"])}')
        p(f'\tSex-Stratified Cohorts: {", ".join(params["sex_strat_cohort_list"])}')
        br()
        p(f'Binary Phenotypes: {", ".join(params["bin_pheno_list"])}')
        p(f'Quantitative Phenotypes: {", ".join(params["quant_pheno_list"])}')
        br()
    
    # Main input files
    with div():
        h2('Main Input Files')
        p(f'Phenotype and Covariate csv: {params["data_csv"]}')
        p(f'Cohort Membership csv: {params["cohort_sets"]}')
        p(f'Sample ID Column: {params["id_col"]}')
    if params["use_sparse_GRM"]:
        with div():
            h2('Sparse GRM Input Files:')
            p(f'Sparse GRM: {params["step1_sparse_grm"]}')
            p(f'Sparse GRM Samples: {params["step1_sparse_grm_samples"]}')
    else:
        with div():
            h2('Pre-Step 1 Plink QC')
            p(f'Plink Prefix: {params["exome_plink_prefix"]}')
            p(f'Min MAF: {params["maf"]}')
            p(f'Max Missingness per Variant: {params["geno"]}')
            p(f'MAX HWE P-Value: {params["hwe"]}')
    br()

    # SAIGE Step 1 configuration options
    with div():
        h2('SAIGE Step 1 Configuration:')
        p(f'Step 1 Script Path: {params["step1_script"]}')
        if params["use_sparse_GRM"]:
            p(f'Plink Prefix: {params["step1_plink_prefix"]}')
        p(f'Categorical Covariates: {", ".join(params["cat_covars"])}')
        p(f'\tSex-Stratified Categorical Covariates: {", ".join(params["sex_strat_cat_covars"])}')
        p(f'Continuous Covariates: {", ".join(params["cont_covars"])}')
        p(f'\tSex-Stratified Continuous Covariates: {", ".join(params["sex_strat_cont_covars"])}')
        br()

    # SAIGE Step 2 configuration options
    with div():
        h2('SAIGE-GENE (Step 2) Configuration:')
        p(f'Step 2 Script Path: {params["step2_script"]}')
        p(f'Exome Plink Prefix: {params["exome_plink_prefix"]}')
        p(f'Chr-Separated Group File Path: {params["group_file_prefix"]}')
        br()
        p(f'Min MAF: {params["min_maf"]}')
        p(f'Min MAC: {params["min_mac"]}')
        p(f'Use Firth?: {params["use_firth"]}')
        p(f'Firth Cutoff: {params["firth_cutoff"]}')
        p(f'Leave-One-Chr-Out (LOCO): {params["LOCO"]}')
        br()
        h3('Gene Burden Variant Criteria:')
        p('Max MAF Cutoffs:')
        with ul():
            # Split these at commas to get a list
            for cutoff in params["grouptest_maf"].split(','):
                li(cutoff)
        p('Annotation Groups:')
        with ul():
            # Split these at commas to get a list
            for group in params["grouptest_annotation"].split(','):
                li(group)
    
    # Post-processing options
    with div():
        h2('Post-Processing Configuration:')
        p(f'Gene Coordinates File: {params["gene_location_file"]}')
        p(f'P Threshold for Top Hits: {params["p_cutoff_summarize"]}')

# Write to HTML file
open('saige_exwas_methods.html', 'w+').write(str(outfile))