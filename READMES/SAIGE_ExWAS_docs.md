
Documentation for SAIGE ExWAS
=============================

# Module Overview


SAIGE ExWAS is a pipeline for doing whole-exome association study of rare variants and gene burdens with traits using SAIGE software

[Paper Link for Reference](https://www.nature.com/articles/s41588-022-01178-w)

[Tool Documentation Link](https://saigegit.github.io/SAIGE-doc/)
## Cloning Github Repository:


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory run commands (our pipelines assume all of the nextflow files/scripts are in the current working directory)
## Software Requirements:


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Command to Pull from Google Container Registry: `docker pull gcr.io/ritchie-aou-psom-9015/saige:latest`

* Run Command: `nextflow run workflows/saige_exwas.nf -profile cluster`

* Common `nextflow run` flags:

    * `-resume` flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning

    * `-stub` performs a sort of dry run of the whole workflow, checks channels without executing any code

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile standard` uses the docker image to executes the processes

    * `-profile cluster` uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems

    * `-profile all_of_us` uses the docker image to execute pipelines on the All of Us Researcher Workbench

* for more information visit the [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Configuration Parameters and Input File Descriptions

## Workflow


* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `quant_pheno_list` (Type: List)

    * Quantitative phenotype list
## Pre-Processing


* `id_col ` (Type: String)

    * ID column label

* `data_csv` (Type: File Path)

    * A csv table with all of the phenotypes and covariates to be tested

    * Corresponding Input File: Phenotypes and Covariates

        * table with participants as rows and all needed phenotypes and covariates as columns

        * Type: Data Table

        * Format: csv

        * Input File Header:





        ```
        IID,y_quantitative,y_binary,x1,x2,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
        1a1,2.0046544617651,0,1.51178116845085,1,0,0,0,0,0,0,0,0,1,0
        1a2,0.104213400269085,0,0.389843236411431,1,0,0,0,0,0,0,0,0,1,1
        1a3,-0.397498354133647,0,-0.621240580541804,1,0,0,0,0,0,0,0,0,0,1
        1a4,-0.333177899030597,0,-2.2146998871775,1,0,0,0,0,0,0,0,0,1,1
        1a5,1.21333962248852,0,1.12493091814311,1,0,0,0,0,0,0,0,0,1,0
        1a6,-0.275411643032321,0,-0.0449336090152309,1,0,0,0,0,0,0,0,0,1,0
        1a7,0.438532936074923,0,-0.0161902630989461,0,0,0,0,0,0,0,0,0,0,0
        1a8,0.0162938047248591,0,0.943836210685299,0,0,0,0,0,0,0,0,0,1,1
        1a9,0.147167262428064,0,0.821221195098089,1,0,0,0,0,0,0,0,0,1,0
        ```

* `cohort_sets` (Type: File Path)

    * A binary csv table in which the columns are the cohorts and the rows are the individuals. A 1 means that individual is a member of the column’s cohort, and a 0 means they aren’t.

    * Corresponding Input File: Cohort Membership

        * 0/1 table with cohorts as columns and participants as rows - 1 indicates that that row’s participant is a member of that column’s cohort

        * Type: Data Table

        * Format: csv

        * Input File Header:





        ```
        IID,POP1,POP2,POP3
        1a1,1,0,1
        1a2,1,0,0
        1a3,0,0,0
        1a4,1,0,0
        1a5,1,0,0
        1a6,1,1,0
        1a7,0,0,0
        1a8,1,0,1
        1a9,0,1,1
        
        ```
## QC Options


* `min_maf` (Type: Float)

    * Minimum minor allele frequency for plink QC
## Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
## SAIGE Step 1


* `geno` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants, genotype rate  filters out all variants with missing call rates exceeding the provided value 

* `step1_script  ` (Type: File Path)

    * Fits the null logistic/linear mixed model using a full or a sparse genetic relationship matrix (GRM). The GRM estimate the genetic relationship between two individuals over a certain number of SNPs

* `exome_plink_prefix` (Type: Plink Fileset Prefix)

    * Exome plink input files 

    * Corresponding Input File: SAIGE Exome Plink Files

        * a hard-call plink set of exome data

        * Type: Plink Set

        * Format: plink binary

        * Input File Header:





        ```
        PMBB-Release-2020-2.0_genetic_exome_GL_norm{.bed,.bim,.fam,.log,.pgen,.psam,.pvar}
        ```

* `group_file_prefix` (Type: Chr File Prefix)

    * Has the variant positions for each gene as well as the variant annotation for each variant in the gene in SAIGE format 

    * Corresponding Input File: SAIGE Group Annotation Files

        * text files formatted like this example from the SAIGE github: 

        * Type: Data Table

        * Format: saige group (txt)

        * Input File Header:





        ```
        ENSG00000000457 var     1_169853716_C_A 1_169853716_C_T 1_169853717_C_CAGTT
        ENSG00000000457 anno    other_missense  damaging_missense       damaging_missense
        ENSG00000000460 var     1_169795119_C_T 1_169795121_G_C 1_169795123_C_G
        ENSG00000000460 anno    other_missense  other_missense  other_missense
        ```

* `hwe` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants

* `step1_sparse_grm_samples` (Type: File Path)

    * List of IDs to use in the sparse GRM 

    * Corresponding Input File: SAIGE Sparse GRM Sample IDs

        * (optional) sample IDs for a sparse relatedness matrix

        * Type: List File

        * Format: txt

        * Input File Header:





        ```
        1a1
        1a2
        1a3
        1a4
        1a5
        1a6
        1a7
        1a8
        1a9
        1a10
        ```

* `maf ` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
## SAIGE Step 2


* `min_mac ` (Type: Float)

    * SAIGE-GENE Step 2 Parameters

* `firth_cutoff` (Type: Float)

    * P-value ()

* `use_firth` (Type: Bool (R: TRUE or FALSE))

    * True to use firth logistic regression

* `grouptest_annotation` (Type: String)

    * Annotations for variants

* `grouptest_maf` (Type: String)

    * MAF cutoffs

* `LOCO` (Type: Bool (R: TRUE or FALSE))

    * Usually a GWAS method 
## Post-Processing


* `region_col_names` (Type: Map (Dictionary))

    * Default SAIGE Region column names mapped to new ones

* `p_cutoff_summarize` (Type: Float)

    * P-Value Threshold for Summarizing Results at the End, arbitrary p-value threshold for creating a table of results combined with low p-values 
## Plotting


* `gene_location_file ` (Type: File Path)

    * This file is used for getting gene-based coordinates for plotting  

    * Corresponding Input File: Gene Location File

        * CSV file of 

        * Type: Data Table

        * Format: tsv

        * Input File Header:





        ```
        gene_id chromosome  seq_region_start    seq_region_end  gene_symbol
        GENE1   1   1   90  GS1
        GENE2   2   91  100 GS2
        ```
# Output Files from SAIGE_ExWAS


* Quantitative Phenotype Violin Plots

    * A violin plot for each quantitative phenotype. One file will be generated for each phenotype and all cohorts will be plotted on each. 

    * Type: Summary Plot

    * Format: png

    * Parallel By: Phenotype

* Singles QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for variants.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

    * Parallel By: Cohort, Phenotype

* Phenotype Summary Table

    * A csv file with summary statistics for the phenotypes. Statistics are computed across all cohorts, phenotypes. For binary phenotypes, it also computes: Total Number, number of cases, Number of Controls, and Prevalence, and summary stats. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    COHORT |PHENO     |N  |Controls|Cases|Prevalence         |mean             |std              |min  |25%  |50%  |75%    |max
    AMR_ALL|AAA       |474|464.0   |10.0 |0.02109704641350211|                 |                 |     |     |     |       |
    AMR_ALL|T2D       |546|425.0   |121.0|0.2216117216117216 |                 |                 |     |     |     |       |
    AMR_ALL|BMI_median|531|        |     |                   |29.95419020715631|7.10618710646257 |15.17|24.84|28.69|33.435 |64.465
    AMR_ALL|LDL_median|310|        |     |                   |93.09516129032258|34.95209123947716|9.0  |70.25|90.0 |112.375|291.0
    ```

* Regions Manhattan Plots

    * A dot plot (manhattan plot) of significant gene regions associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Regions QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for gene regions.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

    * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Singles Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the variant (singles) analysis. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Phenotype

    * Output File Header:





    ```
    phenotype|chromosome|base_pair_location|variant_id        |other_allele|effect_allele|effect_allele_count|effect_allele_frequency|missing_rate|beta      |standard_error|t_statistic|variance|p_value   |p_value_na|is_spa_test|allele_freq_case|allele_freq
    T2Diab   |21        |41801254          |21_41801254_TCTG_T|TCTG        |T            |277                |0.0046611              |0.0         |-0.099231 |0.167775      |-3.52526   |35.5258 |0.5542179 |0.5542179 |False      |0.00426841      |0.00474126
    T2Diab   |21        |41801360          |21_41801360_C_T   |C           |T            |41                 |0.00068991             |0.0         |-0.864441 |0.633121      |-3.98228   |5.38237 |0.08606924|0.08606924|False      |0.000297796     |0.000769948
    T2Diab   |21        |41801603          |21_41801603_C_T   |C           |T            |24                 |0.00040385             |0.0         |0.322923  |0.570593      |0.991852   |3.07148 |0.5714322 |0.5714322 |False      |0.000496327     |0.000384974
    T2Diab   |21        |41801645          |21_41801645_G_A   |G           |A            |58                 |0.000975971            |0.0         |0.0167811 |0.35132       |0.135962   |8.10206 |0.9619027 |0.9619027 |False      |0.00109192      |0.000952304
    ```

* Singles Manhattan Plots

    * A dot plot (manhattan plot) of significant variants associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Cohort, Phenotype

* Singles Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Singles (Variant) Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het,n
    EUR_M,BMI_median,1,1203928,1_1203928_G_A,G,A,13,0.000425393,0.0,6.68476,1.49857,2.97666,0.445291,8.167552e-06,,,,,,,,,,,15280.0
    AFR_M,BMI_median,1,6073749,1_6073749_C_T,C,T,14,0.00182292,0.0,8.90012,1.93388,2.37977,0.267386,4.180559e-06,,,,,,,,,,,3840.0
    AFR_F,LDL_median,1,11134343,1_11134343_G_C,G,C,13,0.00127226,0.0,36.8257,8.24963,0.541106,0.0146937,8.047218e-06,,,,,,,,,,,5109.0
    AFR_ALL,LDL_median,1,11828619,1_11828619_C_T,C,T,12,0.000740466,0.0,45.5826,9.2016,0.538359,0.0118106,7.279194e-07,,,,,,,,,,,8103.0
    ```

* Regions Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the gene (regions) analysis if run. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Phenotype

    * Output File Header:





    ```
    phenotype|gene           |annot            |max_maf|p_value           |p_value_burden    |p_value_skat      |beta_burden        |se_burden         |mac   |mac_case|mac_control|rare_var_count|ultrarare_var_count
    T2Diab   |ENSG00000141956|pLoF             |0.0001 |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.001  |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.01   |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|damaging_missense|0.0001 |0.464219450219203 |0.464219450219203 |0.464219450219203 |-0.0110759683619445|0.0151328276810456|52.0  |7.0     |45.0       |0.0           |41.0
    ```

* Binary Phenotype Bar Plots

    * A bar plot for each binary phenotype. One file will be generated for each phenotype and all cohorts will be plotted. 

    * Type: Summary Plot

    * Format: png

    * Parallel By: Phenotype

* Regions Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Regions Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    cohort,phenotype,gene,annot,max_maf,p_value,p_value_burden,p_value_skat,beta_burden,se_burden,mac,mac_case,mac_control,rare_var_count,ultrarare_var_count
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.0001,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.001,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.01,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_M,LDL_median,ENSG00000006530,pLoF,0.01,8.20696106989754e-06,8.20696106989754e-06,8.20696106989754e-06,3.23068908720413,0.724416423633968,3.0,,,0.0,3.0
    ```
# Example Config File Contents


```

params {
    // DATA FILES
    // ----------
    // all default paths are for PMBB WES
    data_csv = "/path/to/directory/cleaned_test_pheno_covars.csv"
    // data_csv = "/path/to/directory/common_phecodes_covariate_ALL.csv"
    
    // cohort sets
    cohort_sets = "/path/to/directory/Exome_sample_table.csv"
    
    // this is for getting gene-based coordinates for plotting
    gene_location_file = "/path/to/directory/homo_sapiens_111_b38.txt"

    // ID column label
    id_col = "PMBB_ID"
    
    // Full list of cohorts (usually ancestry-stratified and/or sex-stratified)
    cohort_list = [
        "PMBB_AMR_ALL", "PMBB_AMR_F","PMBB_AMR_M",
        "PMBB_AFR_ALL", "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_ALL", "PMBB_EUR_F", "PMBB_EUR_M"
        ]
        
    // smaller list of cohorts for testing
    // cohort_list = [
    //     "PMBB_AMR_ALL", "PMBB_AMR_F","PMBB_AMR_M",
    //     "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M"
    //    ]
    
    // subset of cohorts that are female- or male-only which should exclude sex-based covariates
    sex_strat_cohort_list = [
        "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_F", "PMBB_EUR_M",
        ]

    // smaller list of sex stratified for testing
    // sex_strat_cohort_list = [
    //     "PMBB_AMR_F","PMBB_AMR_M",
    //     "PMBB_EAS_F", "PMBB_EAS_M",
    //    ]
    
    // Full list of chromosomes
    chromosome_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    
    // Small list of chormosomes for testing
    // chromosome_list = ["22"]

    // binary and quantitative phenotype [lists] or path to file of newline-separated lists
    // bin_pheno_list = "/path/to/directory/common_phecodes_list.txt"
    bin_pheno_list = ["T2D", "AAA"]
    quant_pheno_list = ["BMI_median", "LDL_median"]
    // sex_specific pheno file - these will be skipped for _ALL cohorts
    sex_specific_pheno_file = "/path/to/directory/phecode_Sex_specific.txt"

    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Exome_PC1", "Exome_PC2", "Exome_PC3", "Exome_PC4"]
    sex_strat_cat_covars = []
    sex_strat_cont_covars = cont_covars

    // NextFlow, Docker, and Singularity OPTIONS
    // ------------------------------------------
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"

    // default paths assume use of the docker container
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"
    step2_script = "/usr/local/bin/step2_SPAtests.R"

    // gpu paramater either ON or OFF, need to set config to -c nextflow_gpu.config
    GPU = 'OFF'
    
    // Minimum numbers for filtering cohort-phenotype combinations
    min_bin_cases = 100
    min_quant_n = 1000

    // Config parameters for using precomputed sparse GRM:
    // use_sparse_GRM = true
    // step 1 path should be the small subset of markers used to fit the GRM
    // step1_plink_prefix = "/path/to/directory/PMBB_exome_random_autosomal_markers"
    // step1_sparse_grm = "/path/to/directory/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
    // step1_sparse_grm_samples = "/path/to/directory/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

    // Config parameters for using real-time FULL GRM:
    use_sparse_GRM = false

    // Genetic Data Inputs:
    // Without a GRM, the exome plink set is used for step 1 because it needs rare variants
    exome_plink_prefix = "/path/to/directory/PMBB-Release-2020-2.0_genetic_exome_GL_norm"
    group_file_prefix = "/path/to/directory/subset."
    
    // Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
    // Current defaults are recommended by GBMI analysis plan
    maf = 0.01
    geno = 0.01
    hwe = 1E-6
    
    // SAIGE-GENE Step 2 Parameters
    // Current defaults are recommended by BRAVA analysis plan
    min_maf = 0
    min_mac = 0.5
    grouptest_maf = "0.0001,0.001,0.01"
    grouptest_annotation = "pLoF,damaging_missense,other_missense,synonymous,pLoF;damaging_missense,pLoF;damaging_missense;other_missense;synonymous"
    use_firth = "TRUE"
    firth_cutoff = 0.1
    LOCO = "FALSE"

    // P-Value Threshold for Summarizing Results at the End
    p_cutoff_summarize = 0.00001

    // Dictionary (Map) with default SAIGE Region column names mapped to new ones
    regions_col_names = [
        Region: 'gene',
        Group: 'annot',
        max_MAF: 'max_maf',
        Pvalue: 'p_value',
        Pvalue_Burden: 'p_value_burden',
        BETA_Burden: 'beta_burden',
        SE_Burden: 'se_burden',
        Pvalue_SKAT: 'p_value_skat',
        MAC: 'mac',
        MAC_case: 'mac_case',
        MAC_control: 'mac_control',
        Number_rare: 'rare_var_count',
        Number_ultra_rare: 'ultrarare_var_count'
    ]

    // Dictionary (Map) with default SAIGE SingleAssoc column names mapped to new ones
    singles_col_names = [
        CHR: 'chromosome',
        POS: 'base_pair_location',
        MarkerID: 'variant_id',
        Allele1: 'other_allele',
        Allele2: 'effect_allele',
        AC_Allele2: 'effect_allele_count',
        AF_Allele2: 'effect_allele_frequency',
        MissingRate: 'missing_rate',
        BETA: 'beta',
        SE: 'standard_error',
        Tstat: 't_statistic',
        var: 'variance',
        'p.value': 'p_value',
        'p.value.NA': 'p_value_na',
        'Is.SPA': 'is_spa_test',
        AF_case: 'allele_freq_case',
        AF_ctrl: 'allele_freq_ctrl',
        N_case: 'n_case',
        N_ctrl: 'n_ctrl',
        N_case_hom: 'n_case_hom',
        N_case_het: 'n_case_het',
        N_ctrl_hom: 'n_ctrl_hom',
        N_ctrl_het: 'n_ctrl_het',
        N: 'n'
    ]

}

```
# Current Dockerfile for the Container/Image


```docker
# build: for python packages, plink, biofilter, and NEAT plots
FROM continuumio/miniconda3 as build

WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

RUN apt-get update \
    # install tools needed to install plink, biofilter, and NEAT plots
    && apt-get install -y --no-install-recommends git wget unzip libtiff5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install python packages needed for pipeline
    && conda install -y -n base -c conda-forge -c bioconda adjustText apsw sqlite dominate bgenix scipy pandas seaborn matplotlib conda-build numpy \
    && conda clean --all --yes \
    # install plink
    && wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240526.zip \
    && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip \
    && unzip plink2_linux_x86_64_20240526.zip \
    && unzip plink_linux_x86_64_20231211.zip \
    && rm -rf plink2_linux_x86_64_20240526.zip plink_linux_x86_64_20231211.zip \
    # symlink libtiff.so.6 to libtiff.so.5 to overcome error
    && ln -s /opt/conda/lib/libtiff.so.6 /opt/conda/lib/libtiff.so.5 \
    # download NEAT-plots
    && git clone https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    # move manhattan plot scripts to /app
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    # install biofilter
    && wget https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app/ \
    && /opt/conda/bin/python setup.py install \
    # make biofilter script executable
    && chmod a+rx /app/biofilter.py \
    # remove NEAT-plots directory and biofilter executable
    && rm -R NEAT-Plots biofilter.tar.gz

# dev: for SAIGE and SAIGE-dependent packages
FROM rocker/tidyverse:4.1.3 as dev

WORKDIR /tmp

RUN apt-get update && \
    # install tools needed to install SAIGE
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends build-essential cmake libopenblas-base python3-pip r-cran-devtools git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install cget \
    # download SAIGE
    && git clone https://github.com/saigegit/SAIGE.git \
    # install SAIGE-dependent R packages
    && Rscript SAIGE/extdata/install_packages.R \
    # install SAIGE
    && R CMD INSTALL SAIGE \
    # move SAIGE scripts into $PATH
    && mv SAIGE/extdata/step1_fitNULLGLMM.R SAIGE/extdata/step2_SPAtests.R SAIGE/extdata/step3_LDmat.R SAIGE/extdata/createSparseGRM.R /usr/local/bin/ \
    # make SAIGE scripts executable
    && chmod -R a+x /usr/local/bin/ \
    # remove SAIGE directory
    && rm -R SAIGE

# main: file image with only necessary packages and scripts
FROM ubuntu:20.04 as main

WORKDIR /app

# copy conda packages and installed scripts to main image
COPY --from=dev /tmp/ /app/
COPY --from=dev /usr/local/ /usr/local/
COPY --from=build /opt/conda/ /opt/conda/
COPY --from=build /app/ /app/

RUN apt-get update \
    # install R into Ubuntu base image
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install NEAT-plots manhattan plotting python packages
    && /opt/conda/bin/conda develop /app/manhattan-plot/ \
    # move plink executables into $PATH
    && mv plink2 plink /usr/bin

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

USER root

```
# Current `nextflow.config` contents


```
// includeConfig '${launchDir}/configs/saige_exwas.config'
// includeConfig '${launchDir}/configs/saige_gene_phewas.config'
includeConfig '${launchDir}/configs/saige_variant_phewas.config'

profiles {

    non_docker_dev {
        process.executor = 'local'
    }

    standard {
        process.executor = 'local'
        process.container = Enter .sif file
        docker.enabled = true
    }

    cluster {
        process.executor = 'lsf'
        process.queue = Enter the Queue name
        executor {
            queueSize=500
        }
        process.memory = '15GB'
    	process.container = Enter .sif file
        singularity.enabled = true
        singularity.runOptions = '-B /project/,/static/'
    }

    all_of_us {
        process.executor = 'google-lifesciences'
        process.memory = '15GB'
        process.container = Enter .sif file
        google.zone = "us-central1-a"
        google.project = 'verma-pmbb-codeworks-psom-bf87' // change to your project id
        google.lifeSciences.debug = true
        google.lifeSciences.network = "network"
        google.lifeSciences.subnetwork = "subnetwork"
        google.lifeSciences.usePrivateAddress = false
        google.lifeSciences.serviceAccountEmail = 'project-service-account@verma-pmbb-codeworks-psom-bf87.iam.gserviceaccount.com' // change to your service email
        google.lifeSciences.copyImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        google.enableRequesterPaysBuckets = true
        workDir='gs://fc-secure-f3e7d01e-18fa-40ba-bb3e-4d7497ba7d5b/work/' // change to your working directory in your workspace bucket
    }
}

```