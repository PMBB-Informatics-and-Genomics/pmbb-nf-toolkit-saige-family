
Documentation for SAIGE-GATE
============================

# Module Overview


nan
- [Tool Paper Link for Reference](nan)
- [Tool Documentation Link for Reference](nan)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gate.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 24.04.3](https://www.nextflow.io/docs/latest/cli.html)

  > **Warning:** Nextflow 25.x and 26.x introduce breaking changes incompatible with this pipeline. Use version 24.04.3 exactly.

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Pull from Google Container Registry: `docker pull gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest`

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_gwas.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up workflow where it left off

    * `-stub` performs a dry run, checks channels without executing code

    * `-profile` selects the compute profiles in nextflow.config

    * `-profile standard` uses the Docker image to execute processes

    * `-profile cluster` uses the Singularity container and submits processes to a queue

    * `-profile all_of_us` uses the Docker image on All of Us Workbench

* More info: [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Detailed Pipeline Steps

## Part I: Setup


1. Start your own tools directory and go there. You may do this in your project analysis directory, but it often makes sense to clone into a general `tools` location

```sh
# Make a directory to clone the pipeline into
TOOLS_DIR="/path/to/tools/directory"
mkdir $TOOLS_DIR
cd $TOOLS_DIR
```

2. Download the source code by cloning from git

```sh
git clone None
cd $TOOLS_DIR/pmbb-nf-toolkit-saige-family
```

3. Build the singularity image
    - you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
    - this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly.


```sh
cd $TOOLS_DIR/pmbb-nf-toolkit-saige-family
singularity build saige.sif docker://pennbiobank/saige:latest
```
## Part II: Configure your run


1. Make a separate analysis/run/working directory.
    - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses.
    - ❗This step is optional, but We Highly recommend making a `tools` directory separate from your `run` directory. We recommend storing the `nextflow.conf` in here as it shouldn't change between runs.


```sh
WDIR="/path/to/analysis/run1"
mkdir -p $WDIR
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/blob/main/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile $MY_PROFILE`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and passing it to the pipeline with `-c /path/to/nextflow.config`.


3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
    - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters.
    - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
    - you can compartamentalize your config file as much as you like by passing
    - There are 2 ways to specify the config file during a run:

        - with the `-c` option on the command line: `nextflow run -c SAIGE_FAMILY/configs/saige_gate.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig SAIGE_FAMILY/configs/saige_gate.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-saige-family/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gate.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gate.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gate.config \
   -resume
```
# Pipeline Parameters

---

## Required Parameters

Every run must set all of the following. Everything else has a default or is optional.

| Parameter | Section | Description |
|-----------|---------|-------------|
| `data_csv` | Input Files | CSV of participant phenotypes and covariates |
| `cohort_sets` | Input Files | Cohort membership table (0/1 by participant) |
| `id_col` | Step 1 | Participant ID column name |
| `cohort_list` | Step 1 | List of cohorts to run |
| `sex_strat_cohort_list` | Step 1 | Sex-stratified subset of `cohort_list` |
| `bin_pheno_list` and/or `quant_pheno_list` and/or `survival_pheno_list` | Step 1 | Phenotype list(s); at least one must be non-empty |
| `cat_covars` | Step 1 | Categorical covariate column names |
| `cont_covars` | Step 1 | Continuous covariate column names |
| `chromosome_list` | Step 1 | Chromosomes to include in Step 2 |
| `step1_plink_prefix` | Input Files | PLINK fileset prefix for Step 1 GRM |
| `ftype` | Input Files | Step 2 genotype format: `”PLINK1”`, `”PLINK2”`, or `”BGEN”` |
| `step2_plink_prefix` / `step2_pgen_prefix` / `step2_bgen_prefix` | Input Files | Step 2 genotype file prefix (matches `ftype`) |

---

## Input Files

These are the data files the pipeline reads. All paths must be accessible from the compute node.

### Phenotype and Cohort Data

| Parameter | Required | Description |
|-----------|----------|-------------|
| `data_csv` | **Required** | Path to a CSV table with one row per participant. Columns are phenotypes and covariates. Must contain the column named by `id_col`. For survival analyses, must also contain the event-time column specified by `event_time_col`. |
| `cohort_sets` | **Required** | Path to a CSV table where columns are cohort names and rows are participants. A `1` means the participant is in that cohort, `0` means they are not. |
| `sex_specific_pheno_file` | Optional | Path to a newline-separated text file listing phenotypes that should only be included in sex-stratified cohorts. Leave as `null` if not used. |

**Example `data_csv` format:**
```
IID,y_quantitative,y_binary,age_at_event,event_occurred,x1,x2
1a1,2.004,0,65.2,1,1.511,1
1a2,0.104,0,72.1,0,0.389,1
```

### Step 1 Genetic Input Files

| Parameter | Required | Description |
|-----------|----------|-------------|
| `step1_plink_prefix` | **Required** | File prefix for the PLINK binary fileset (`.bed/.bim/.fam`) used to build the GRM in Step 1. Should cover all chromosomes in a single merged file. |
| `step1_sparse_grm` | Optional | Path to a precomputed sparse GRM file (`.sparseGRM.mtx`). Only used when `use_sparse_GRM = true`. |
| `step1_sparse_grm_samples` | Optional | Path to a sample ID file for the sparse GRM. Only used when `use_sparse_GRM = true`. |

### Step 2 Genetic Input Files

Step 2 runs the genome-wide association test. Set `ftype` to choose the file format.

| Parameter | Required when | Description |
|-----------|---------------|-------------|
| `ftype` | **Always** | Genotype file format for Step 2. Must be `”PLINK1”`, `”PLINK2”`, or `”BGEN”`. |
| `step2_plink_prefix` | `ftype = “PLINK1”` | Prefix for chromosome-separated PLINK1 files. |
| `step2_pgen_prefix` | `ftype = “PLINK2”` | Prefix for chromosome-separated PLINK2 files. |
| `step2_bgen_prefix` | `ftype = “BGEN”` | Prefix for chromosome-separated BGEN files. |
| `bgen_samplefile` | `ftype = “BGEN”` | Path to the `.sample` file for the BGEN dataset. |

---

## Step 1 — Cohort Setup & Null Model Parameters

Step 1 sets up cohort/phenotype combinations and fits a null GLMM (or Cox model for survival). It also runs PLINK2 QC to select high-quality variants for GRM construction.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `id_col` | *(required)* | The name of the participant ID column in `data_csv` and `cohort_sets`. |
| `cohort_list` | *(required)* | List of cohort names to run. Must match column names in `cohort_sets`. |
| `sex_strat_cohort_list` | *(required)* | Subset of `cohort_list` that are sex-stratified. |
| `bin_pheno_list` | `[]` | List of binary phenotype column names, OR a path to a file with one phenotype per line. |
| `quant_pheno_list` | `[]` | List of quantitative phenotype column names, OR a path to a file with one phenotype per line. |
| `survival_pheno_list` | `[]` | List of survival (time-to-event) phenotype column names, OR a path to a file with one phenotype per line. Leave as `[]` if not running survival analyses. |
| `chromosome_list` | *(required)* | List of chromosomes to include in Step 2. For a full GWAS: `[“1”,”2”,...,”22”]`. |
| `cat_covars` | `[]` | Categorical covariate column names (e.g. `[“SEX”]`). Used for non-sex-stratified cohorts. |
| `cont_covars` | `[]` | Continuous covariate column names (e.g. `[“AGE”,”PC1”,”PC2”]`). Used for non-sex-stratified cohorts. |
| `use_sparse_GRM` | `false` | Set to `true` to use a precomputed sparse GRM. Requires `step1_sparse_grm` and `step1_sparse_grm_samples`. |
| `maf` | `0.01` | Minimum minor allele frequency for Step 1 variant QC. |
| `geno` | `0.01` | Maximum per-variant missingness for Step 1 QC. |
| `hwe` | `1E-6` | Hardy-Weinberg equilibrium p-value threshold for Step 1 QC. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sex_strat_cat_covars` | `[]` | Categorical covariates for sex-stratified cohorts. Use to drop sex from covariates for male-only or female-only cohorts. |
| `sex_strat_cont_covars` | same as `cont_covars` | Continuous covariates for sex-stratified cohorts. |
| `min_bin_cases` | `50` | Minimum number of cases required to run a binary phenotype for a given cohort. |
| `min_quant_n` | `500` | Minimum number of non-missing samples required to run a quantitative phenotype for a given cohort. |
| `min_survival_cases` | `100` | Minimum number of events required to run a survival phenotype for a given cohort. |
| `event_time_col` | *(required for survival)* | Name of the event-time column in `data_csv`. Only needed when `survival_pheno_list` is non-empty. |
| `event_time_bin` | *(required for survival)* | Bin size for discretizing event times in SAIGE's survival model. Use the same unit as `event_time_col` (e.g. `1` for years, `1/12` for months). |
| `LOCO` | `”TRUE”` | Leave-one-chromosome-out for the null model. Recommended `”TRUE”` for GWAS. |
| `GPU` | `”OFF”` | Set to `”ON”` to use GPU-accelerated SAIGE-DOE for Step 1. |
| `step1_script` | `/usr/local/bin/step1_fitNULLGLMM.R` | Path to the SAIGE Step 1 R script. Only change if using a custom container. |

---

## Step 2 — Association Test Parameters

Step 2 tests each variant genome-wide for association using the null model from Step 1.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_mac` | `40` | Minimum minor allele count for a variant to be tested. |
| `min_maf` | `0` | Minimum minor allele frequency for a variant to be tested. Typically left at `0` and controlled by `min_mac`. |
| `use_firth` | `TRUE` | Whether to use Firth logistic regression as a fallback for small p-values in binary traits. |
| `firth_cutoff` | `0.1` | P-value threshold below which Firth regression is applied. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LOCO` | `”TRUE”` | Leave-one-chromosome-out for Step 2. Should match the Step 1 setting. |
| `isImputed` | `”FALSE”` | Set to `”TRUE”` if Step 2 genotype files are imputed (e.g. BGEN from TOPMed). |
| `minInfo` | `0` | Minimum imputation info score. Only active when `isImputed = “TRUE”`. Typical value: `0.3`. |
| `step2_script` | `/usr/local/bin/step2_SPAtests.R` | Path to the SAIGE Step 2 R script. Only change if using a custom container. |

---

## Post-Processing Parameters

These parameters control result filtering, annotation, and visualization.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `annotate` | `true` | Set to `true` to annotate results with RSIDs and nearest gene information using Biofilter. |
| `p_cutoff_summarize` | `0.00001` | P-value threshold for the top-hits summary table. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gene_location_file` | `/app/NCBI.gene.loc` | Path to a gene coordinate file for Manhattan plot labeling. Bundled in the container. |
| `biofilter_loki` | *(required if annotating)* | Path to a `loki.db` database file for Biofilter gene annotations. |
| `biofilter_build` | *(required if annotating)* | Genome build for Biofilter: `”19”` for GRCh37 or `”38”` for GRCh38. |
| `biofilter_close_dist` | `5E4` | Distance in base pairs defining “close” vs. “far” for nearest-gene annotation. Default is 50 kb. |
| `gwas_col_names` | *(see example config)* | Map renaming SAIGE's default output column headers to preferred names. |

---

## Infrastructure (Advanced)

These parameters set executable paths and environment-specific behavior. Defaults are correct for the SAIGE container.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `my_python` | `/opt/conda/bin/python` | Path to the Python executable. |
| `my_bgenix` | `/opt/conda/bin/bgenix` | Path to the `bgenix` executable. Only needed when `ftype = “BGEN”`. |
| `host` | `””` | Execution environment hint. Set to `”AOU”`, `”DNAnexus”`, `”LPC”`, or leave empty. |

---

# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"
    my_bgenix = "/opt/conda/bin/bgenix"

    //setting file type for step 2 (PLINK/BGEN)
    //ftype = "PLINK"
    ftype = "BGEN"
    GPU="OFF"
    annotate=true   
    
    use_sparse_GRM = false
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"
    step2_script = "/usr/local/bin/step2_SPAtests.R"

    data_csv = "/path/to/data/common_ICD_covariate_ALL.csv"

    cohort_sets = "/path/to/data/Imputed_sample_table.csv"

    // default paths are for BIOBANK Geno data
    step1_plink_prefix  = "/path/to/data/pruned_data"
    step2_plink_prefix = "/path/to/data/pruned_data"
    
    // default paths for Imputed Geno data BGEN
    step2_bgen_prefix  = "/path/to/data/BIOBANK-Release-VERSION_genetic_imputed-topmed-r2_chr"
    bgen_samplefile = "/path/to/data/BIOBANK-Release-VERSION_genetic_imputed-topmed-r2_bgen.sample"
    
    
    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Genotype_PC1","Genotype_PC2","Genotype_PC3",
                   "Genotype_PC4", "Genotype_PC5","Genotype_PC6","Genotype_PC7",
                   "Genotype_PC8","Genotype_PC9","Genotype_PC10"]

    min_bin_cases = null
    min_quant_n = null

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
    
    thin_count = ""
    host = "LPC"

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
    bin_pheno_list = "/path/to/data/common_ICD_list.txt"
    quant_pheno_list = []

    // use only for time-to-event gwas 
    survival_pheno_list = []
    event_time_col = "" 
    event_time_bin =  // event time bin size 
    // (must be in same unit of measure of event time col (months,years,days))
    // (e.g. event_time_col = ???? )
    sex_specific_pheno_file = "/path/to/data/icd_Sex_specific.txt"
    
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
## Current `nextflow.config` contents


```
//includeConfig 'configs/saige_exwas.config'
//includeConfig 'configs/saige_gene_phewas.config'
//includeConfig 'configs/saige_variant_phewas.config'

profiles {
    non_docker_dev {
        // run locally without docker
        process.executor = awsbatch-or-lsf-or-slurm-etc
    }

    standard {
        // run locally with docker
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.container = 'pennbiobank/saige:latest'
        docker.enabled = true
    }

    cluster {
        // run on LSF cluster
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        executor {
            queueSize=500
        }
        process.memory = '15GB'
    	process.container = 'saige.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }

    all_of_us {
        // CHANGE EVERY TIME! These are specific for each user, see docs
        google.lifeSciences.serviceAccountEmail = service@email.gservicaaccount.com
        workDir = /path/to/workdir/ // can be gs://
        google.project = terra project id

        // These should not be changed unless you are an advanced user
        process.container = 'gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest' // GCR SAIGE docker container (static)

        // these are AoU, GCR parameters that should NOT be changed
        process.memory = '15GB' // minimum memory per process (static)
        process.executor = awsbatch-or-lsf-or-slurm-etc
        google.zone = "us-central1-a" // AoU uses central time zone (static)
        google.location = "us-central1"
        google.lifeSciences.debug = true 
        google.lifeSciences.network = "network"
        google.lifeSciences.subnetwork = "subnetwork"
        google.lifeSciences.usePrivateAddress = false
        google.lifeSciences.copyImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
        google.enableRequesterPaysBuckets = true
        // google.lifeSciences.bootDiskSize = "20.GB" // probably don't need this
        process{
                withName: 'call_saige_step1_bin' {
                container = '/tnnandi/saige-doe:2'
                }
                }
        process{
                withName: 'call_saige_step1_quant' {
                container = '/tnnandi/saige-doe:2'
                }
        }

    }   
    dnanexus{
        //This is a profile for running on a single machine on dnanexus
        process {
            container = 'pennbiobank/saige:latest'
            executor = 'local'
            memory = '15GB'
        }
        docker {
            enabled = true
        }
        process{
            withName: 'call_saige_step1_bin' {
            container = '/tnnandi/saige-doe:2'
            }
        }
        process{
            withName: 'call_saige_step1_quant' {
            container = '/tnnandi/saige-doe:2'
            }
        }
    }
}

params {
    skip_postprocessing_errors = false
}

process {
    withLabel: safe_to_skip {
        errorStrategy=params.skip_postprocessing_errors ? 'ignore' : 'terminate'
    }
}

```