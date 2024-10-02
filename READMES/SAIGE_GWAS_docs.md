
Documentation for SAIGE GWAS
============================

# Module Overview


SAIGE GWAS is a pipeline for performing genome wide association studies of variants using the R-based SAIGE software. This module has the option of using the biofilter database to provide nearest gene annotation.

[Paper Link for Reference](https://www.nature.com/articles/s41588-018-0184-y)

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

* Run Command: `nextflow run workflows/saige_gwas.nf -config nextflow.config -profile cluster`

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


* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `quant_pheno_list` (Type: List)

    * Quantitative phenotype list

* `cohort_list` (Type: List)

    * List of cohorts usually ancestry stratified and or sex stratified

* `chromosome_list` (Type: List)

    * List of chromosomes, for testing use smaller chromosomes e.g chromosome_list = ["20", "21", "22"]

* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List OR File Path)

    * file path to list of binary phenotypes
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

* `sex_specific_pheno_file` (Type: File Path)

    * A newline-separated list of phenotypes that should only be included in sex-stratified cohorts (e.g., AFR_F but not AFR_ALL)

    * Corresponding Input File: Sex Specific Phenotype List

        * A newline-separated list of phenotypes that should be excluded from non-sex-stratified cohorts (e.g., AFR_F or AFR_M but not AFR_ALL)

        * Type: List File

        * Format: txt
## Association Test Modeling


* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges

* `sex_strat_cont_covars` (Type: List)

    * Continuous covariates for sex stratified cohorts to ensure model converges
## SAIGE Step 1


* `step1_script  ` (Type: File Path)

    * Fits the null logistic/linear mixed model using a full or a sparse genetic relationship matrix (GRM). The GRM estimate the genetic relationship between two individuals over a certain number of SNPs

* `step1_plink_prefix` (Type: Plink Fileset Prefix)

    * Step1 exome plink input fileset  - should be all chromosomes together

    * Corresponding Input File: SAIGE Step 1 Plink Files

        * a hard-call plink set to use for step 1 (usually also exome or genotype files)

        * Type: Plink Set

        * Format: plink binary

        * Input File Header:





        ```
        nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{bed,bim,fam}
        ```

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

* `step1_sparse_grm` (Type: File Path)

    * File Path to precomputed Sparse GRM

    * Corresponding Input File: SAIGE Sparse GRM

        * (optional) a sparse relatedness matrix for faster step 1 computation

        * Type: GRM

        * Format: R sparse matrix

        * Input File Header:





        ```
        sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
        ```
## SAIGE Step 2


* `ftype` (Type: String)

    * “PLINK” or “BGEN” based on input type for steps 1 & 2 

* `step2_script` (Type: File Path)

    * Performs region or gene-based association tests    

* `use_firth` (Type: Bool (R: TRUE or FALSE))

    * True to use firth logistic regression
## Post-Processing


* `p_cutoff_summarize` (Type: Float)

    * P-Value Threshold for Summarizing Results at the End, arbitrary p-value threshold for creating a table of results combined with low p-values 

* `annotate` (Type: Bool (Java: true or false))

    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.

* `biofilter_build` (Type: String)

    * The build to pass to biofilter - can be 19 or 38

* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4

* `biofilter_loki` (Type: File Path)

    * The path to a loki.db file to be used for nearest gene annotation

* `biofilter_script` (Type: File Path)

    * The path to the biofilter script to use. If using the singularity container, should be ‘/app/biofilter.py’
# Output Files from SAIGE_GWAS

# Example Config File Contents


```
params {
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"
    my_bgenix = "/opt/conda/bin/bgenix"

    //setting file type for step 2 (PLINK/BGEN)
    //ftype = "PLINK"
    ftype = "BGEN"
    GPU="ON"
    annotate=true   

    use_sparse_GRM = false
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"
    step2_script = "/usr/local/bin/step2_SPAtests.R"

    data_csv = "/path/to/directory/common_ICD_covariate_ALL.csv"

    cohort_sets = "/path/to/directory/Imputed_sample_table.csv"

    // default paths are for PMBB Geno data
    step1_plink_prefix  = "/path/to/directory/pruned_data"
    step2_plink_prefix = "/path/to/directory/pruned_data"
    
    // default paths for Imputed Geno data BGEN
    step2_bgen_prefix  = "/path/to/directory/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr"
    bgen_samplefile = "/path/to/directory/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_bgen.sample"
    
    
    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Genotype_PC1","Genotype_PC2","Genotype_PC3",
                   "Genotype_PC4", "Genotype_PC5","Genotype_PC6","Genotype_PC7",
                   "Genotype_PC8","Genotype_PC9","Genotype_PC10"]


    sex_strat_cat_covars = []
    sex_strat_cont_covars = cont_covars

    // P-Value Threshold for Summarizing Results at the End
    p_cutoff_summarize = 0.00001

    // ID column label
    id_col = "PMBB_ID"

    // Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
    // Current defaults are recommended by GBMI analysis plan
    maf = 0.01
    geno = 0.01
    hwe = 1E-6

   //Step 2 Parameters
    min_maf = 0
    min_mac = 40
    firth_cutoff = 0.1
    LOCO = "TRUE"
    inverseNormalize="TRUE"

 // this is for getting gene-based coordinates for plotting
    // also wrapped in the docker container
    gene_location_file = "/app/NCBI.gene.loc"
         
    cohort_list = [
        "PMBB_AMR_ALL", "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_ALL", "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_ALL", "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_ALL", "PMBB_SAS_F", "PMBB_SAS_M",
        ]

    // subset of cohorts that are female- or male-only which should exclude sex-based covariates
    sex_strat_cohort_list = [
        "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_F", "PMBB_SAS_M"
        ]

    // binary and quantitative phenotype lists
    bin_pheno_list = "/path/to/directory/common_ICD_list.txt"
    quant_pheno_list = []

    sex_specific_pheno_file = "/path/to/directory/icd_Sex_specific.txt"
    
    gwas_col_names = [
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
        N_ctrl_het: 'n_ctrl_het'
    ]

    // list of chromosomes
     chromosome_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
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
# Advanced Nextflow Users: Take/Emit Info

## Input Channel (take) Description


NONE
## Output Channel (emit) Description


Singles_merge_output- a Channel of two paths to the merged raw and filtered output from SAIGE Step 2
Pheno_table- a Channel for a file that contains sample informations on each phenotype