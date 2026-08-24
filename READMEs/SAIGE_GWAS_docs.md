
Documentation for SAIGE GWAS
============================

# Module Overview


SAIGE GWAS is a pipeline for performing genome wide association studies of variants using the R-based SAIGE software. This module has the option of using the biofilter database to provide nearest gene annotation. 

Please see 
- [Tool Paper Link for Reference](https://www.nature.com/articles/s41588-018-0184-y)
- [Tool Documentation Link for Reference](https://saigegit.github.io/SAIGE-doc/)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gwas.config)
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

        - with the `-c` option on the command line: `nextflow run -c SAIGE_FAMILY/configs/saige_gwas.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig SAIGE_FAMILY/configs/saige_gwas.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-saige-family/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config \
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
| `id_col` | Pre-Processing | Participant ID column name in `data_csv` and `cohort_sets` |
| `cohort_list` | Pre-Processing | List of cohorts to run |
| `sex_strat_cohort_list` | Pre-Processing | Sex-stratified subset of `cohort_list` |
| `bin_pheno_list` and/or `quant_pheno_list` | Pre-Processing | Phenotype list(s); at least one must be non-empty |
| `cat_covars` | Pre-Processing | Categorical covariate column names |
| `cont_covars` | Pre-Processing | Continuous covariate column names |
| `chromosome_list` | Pre-Processing | Chromosomes to include in Step 2 |
| `step1_plink_prefix` | Input Files | PLINK fileset prefix for Step 1 GRM construction |
| `ftype` | Input Files | Step 2 genotype format: `"PLINK1"`, `"PLINK2"`, or `"BGEN"` |
| `step2_plink_prefix` / `step2_pgen_prefix` / `step2_bgen_prefix` | Input Files | Step 2 genotype file prefix (matches `ftype`) |
| `biofilter_loki` | Post-Processing | Path to `loki.db` (required when `annotate = true`) |
| `biofilter_build` | Post-Processing | Genome build for Biofilter: `"19"` or `"38"` (required when `annotate = true`) |

---

## Input Files

These are the data files the pipeline reads. All paths must be accessible from the compute node.

### Phenotype and Cohort Data

| Parameter | Required | Description |
|-----------|----------|-------------|
| `data_csv` | **Required** | Path to a CSV table with one row per participant. Columns are phenotypes and covariates. Must contain the column named by `id_col`. |
| `cohort_sets` | **Required** | Path to a CSV table where columns are cohort names and rows are participants. A `1` means the participant is a member of that cohort, `0` means they are not. The first column must be the participant ID matching `id_col`. |
| `sex_specific_pheno_file` | Optional | Path to a newline-separated text file listing phenotypes that should only be included in sex-stratified cohorts (e.g. `EUR_F`, `EUR_M`) and excluded from combined cohorts (e.g. `EUR_ALL`). Leave as `null` if not used. |

**Example `data_csv` format:**
```
IID,y_quantitative,y_binary,x1,x2,a1,a2
1a1,2.004,0,1.511,1,0,0
1a2,0.104,0,0.389,1,0,0
```

**Example `cohort_sets` format:**
```
IID,POP1,POP2,POP3
1a1,1,0,1
1a2,1,0,0
```

### Step 1 Genetic Input Files

Step 1 fits the null model. It requires a set of high-quality, approximately LD-independent SNPs.

| Parameter | Required | Description |
|-----------|----------|-------------|
| `step1_plink_prefix` | **Required** | File prefix for the PLINK binary fileset (`.bed/.bim/.fam`) used to build the GRM in Step 1. Should cover all chromosomes in a single merged file. The pipeline appends `.bed`, `.bim`, `.fam` automatically. |
| `step1_sparse_grm` | Optional | Path to a precomputed sparse GRM file (`.sparseGRM.mtx`). Only used when `use_sparse_GRM = true`. |
| `step1_sparse_grm_samples` | Optional | Path to a text file with one sample ID per line, listing the samples in `step1_sparse_grm`. Only used when `use_sparse_GRM = true`. |

**Example `step1_plink_prefix`:**
```
step1_plink_prefix = "/path/to/data/pruned_data"
# Pipeline will look for: /path/to/data/pruned_data.bed, .bim, .fam
```

### Step 2 Genetic Input Files

Step 2 runs the association test. It requires chromosome-separated genetic files. Set `ftype` to choose the format.

| Parameter | Required when | Description |
|-----------|---------------|-------------|
| `ftype` | **Always** | Genotype file format for Step 2 (and used to determine which plink set to use for Step 1 sample QC). Must be one of: `"PLINK1"`, `"PLINK2"`, or `"BGEN"`. Case-sensitive. |
| `step2_plink_prefix` | `ftype = "PLINK1"` | Prefix for chromosome-separated PLINK1 files. The chromosome number is appended directly, e.g. prefix `"data.chr"` → `"data.chr1.bed"`. |
| `step2_pgen_prefix` | `ftype = "PLINK2"` | Prefix for chromosome-separated PLINK2 files (`.pgen/.pvar/.psam`). The chromosome number is appended directly. |
| `step2_bgen_prefix` | `ftype = "BGEN"` | Prefix for chromosome-separated BGEN files (`.bgen` + `.bgen.bgi` index). The chromosome number is appended directly. |
| `bgen_samplefile` | `ftype = "BGEN"` | Path to the `.sample` file associated with the BGEN dataset, used to map participant IDs. |

**Example BGEN prefix:**
```
step2_bgen_prefix = "/path/to/data/imputed_chr"
# Pipeline will look for: /path/to/data/imputed_chr1.bgen, imputed_chr1.bgen.bgi, etc.
```

---

## Cohort & Phenotype Setup

These parameters define which cohorts, phenotypes, and covariates are used across both Step 1 and Step 2.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `id_col` | *(required)* | The name of the participant ID column in `data_csv` and `cohort_sets`. E.g. `"IID"` or `"PMBB_ID"`. |
| `cohort_list` | *(required)* | List of cohort names to run. These must match column names in `cohort_sets`. E.g. `["EUR_ALL", "EUR_F", "EUR_M"]`. |
| `sex_strat_cohort_list` | *(required)* | Subset of `cohort_list` that are sex-stratified (female-only or male-only). These cohorts will use `sex_strat_cat_covars` and `sex_strat_cont_covars` instead of the default covariate lists. E.g. `["EUR_F", "EUR_M"]`. |
| `bin_pheno_list` | `[]` | List of binary phenotype column names from `data_csv`, OR a path to a text file with one phenotype per line. Leave as `[]` if running only quantitative phenotypes. |
| `quant_pheno_list` | `[]` | List of quantitative phenotype column names from `data_csv`, OR a path to a text file with one phenotype per line. Leave as `[]` if running only binary phenotypes. |
| `survival_pheno_list` | `[]` | List of survival phenotype column names from `data_csv`, OR a path to a text file with one phenotype per line. Leave as `[]` if not running survival analyses. |
| `chromosome_list` | *(required)* | List of chromosomes to include in Step 2. For a full GWAS: `["1","2",...,"22"]`. For testing, use a small subset like `["20","21","22"]`. |
| `cat_covars` | `[]` | List of categorical covariate column names (e.g. `["SEX"]`). Used for all cohorts not in `sex_strat_cohort_list`. |
| `cont_covars` | `[]` | List of continuous covariate column names (e.g. `["AGE","PC1","PC2"]`). Used for all cohorts not in `sex_strat_cohort_list`. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sex_strat_cat_covars` | `[]` | Categorical covariates for sex-stratified cohorts only. Use this to drop sex from covariates for male-only or female-only cohorts. |
| `sex_strat_cont_covars` | same as `cont_covars` | Continuous covariates for sex-stratified cohorts. Typically the same as `cont_covars`. |
| `min_bin_cases` | `50` | Minimum number of cases required to run a binary phenotype for a given cohort. Cohort/phenotype combinations with fewer cases are skipped. Set to `null` to use the default of 50. |
| `min_quant_n` | `500` | Minimum number of non-missing samples required to run a quantitative phenotype for a given cohort. Cohort/phenotype combinations with fewer samples are skipped. Set to `null` to use the default of 500. |
| `min_survival_cases` | `100` | Minimum number of events required to run a survival phenotype for a given cohort. |

---

## SAIGE Step 1 Parameters

Step 1 fits a null generalized linear mixed model (GLMM) per cohort/phenotype combination. It also runs a PLINK2 QC step to select high-quality variants for GRM construction.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `use_sparse_GRM` | `false` | Set to `true` to use a precomputed sparse GRM (faster for large cohorts). Requires `step1_sparse_grm` and `step1_sparse_grm_samples` to be set. Set to `false` to compute the GRM from the Step 1 plink data. |
| `LOCO` | `"TRUE"` | Leave-one-chromosome-out (LOCO) for the null model. Recommended to keep as `"TRUE"` for GWAS to avoid proximal contamination. Use `"FALSE"` for testing or ExWAS. |
| `maf` | `0.01` | Minimum minor allele frequency for Step 1 variant QC (plink2 `--maf`). Variants below this threshold are excluded from GRM construction. |
| `geno` | `0.01` | Maximum per-variant missingness for Step 1 variant QC (plink2 `--geno`). Variants above this threshold are excluded. |
| `hwe` | `1E-6` | Hardy-Weinberg equilibrium p-value threshold for Step 1 QC (plink2 `--hwe`). Variants with HWE p-value below this threshold are excluded. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_vars_for_GRM` | `null` (cap at 150k) | Maximum number of variants to use for GRM construction. If `null`, the pipeline will error if more than 150,000 variants pass QC, forcing you to either set this parameter or reduce the marker set. If set, the pipeline randomly subsamples down to this number. |
| `min_vars_for_GRM` | `30000` | Minimum number of variants required after QC to proceed with GRM construction. The pipeline will error if fewer variants remain. |
| `min_rare_vars_for_GRM` | `300` | (SAIGE-GENE only) Minimum number of rare variants (MAC 10-20 and MAC >20) required for variance ratio estimation. Not typically relevant for GWAS. |
| `pruning_r2_for_GRM` | `0.6` | LD pruning r² threshold used in the plink2 `--indep-pairwise` step before GRM construction. |
| `GPU` | `"OFF"` | Set to `"ON"` to use GPU-accelerated SAIGE-DOE for Step 1. Requires special hardware and container. Leave as `"OFF"` for standard runs. |
| `step1_script` | `/usr/local/bin/step1_fitNULLGLMM.R` | Path to the SAIGE Step 1 R script inside the container. Only change if using a custom container. |
| `event_time_col` | *(required for survival)* | Name of the event-time column in `data_csv` for survival analyses. Only needed when `survival_pheno_list` is non-empty. |
| `event_time_bin` | *(required for survival)* | Bin size for discretizing event times in SAIGE's survival model. Only needed when `survival_pheno_list` is non-empty. |

---

## SAIGE Step 2 Parameters

Step 2 tests each variant genome-wide for association using the null model from Step 1.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_mac` | `40` | Minimum minor allele count (MAC) for a variant to be tested. Variants below this threshold are skipped. Recommended: 20 for array data, 40 for imputed data. |
| `min_maf` | `0` | Minimum minor allele frequency for a variant to be tested. Typically left at `0` and controlled by `min_mac` instead. |
| `use_firth` | `TRUE` | Whether to use Firth logistic regression as a fallback for variants with small p-values in binary traits. Strongly recommended to keep as `TRUE`. |
| `firth_cutoff` | `0.1` | P-value threshold below which Firth regression is applied in Step 2 (BGEN runs). Only applies to binary phenotypes. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LOCO` | `"TRUE"` | Leave-one-chromosome-out for Step 2. Should match the setting used in Step 1. |
| `isImputed` | `"FALSE"` | Set to `"TRUE"` if your Step 2 genotype files are imputed (e.g. BGEN from TOPMed). When `"TRUE"`, SAIGE uses the imputation info score as a filter. |
| `minInfo` | `0` | Minimum imputation info score for variants to be tested. Only active when `isImputed = "TRUE"`. Typical value: `0.3`. |
| `step2_script` | `/usr/local/bin/step2_SPAtests.R` | Path to the SAIGE Step 2 R script inside the container. Only change if using a custom container. |
| `enable_chunking` | `false` | Set to `true` to split each chromosome into sub-chromosome chunks for Step 2. Useful for very large datasets or fine-grained parallelism. Requires `chunks_manifest` to be set. |
| `chunks_manifest` | *(required if chunking)* | Path to a tab-separated manifest file mapping each chunk to its chromosome and full file path stem (no extension). Format: `chromosome<TAB>/full/path/to/file_stem`. One row per chunk. Chromosome values must match entries in `chromosome_list`. Paths may be absolute or relative to the manifest file's directory. Generate with: `ls *.pgen \| sed 's/.pgen$//' \| awk '{match($0, /_chr([0-9XY]+)_/, m); print m[1] "\t" $0}' > manifest.tsv` |
| `gwas_col_names` | *(see example config)* | Map renaming SAIGE's default output column headers to your preferred names. Only change if downstream tools expect specific column names. |

---

## Post-Processing and Annotation Parameters

These parameters control result filtering, annotation with nearest genes, and visualization.

### Beginner (Most Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `annotate` | `true` | Set to `true` to annotate results with RSIDs and nearest gene information using Biofilter. Set to `false` to skip annotation and produce unannotated Manhattan/QQ plots. |
| `p_cutoff_summarize` | `0.00001` | P-value threshold for the summary table of top hits. Variants with p-value below this threshold are included in the combined summary output file. |

### Advanced (Less Commonly Changed)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `biofilter_script` | `/app/biofilter.py` | Path to the Biofilter Python script. The default path is correct when using the SAIGE container. Only change if running Biofilter outside the container. |
| `biofilter_loki` | *(required if annotating)* | Path to a `loki.db` database file used by Biofilter for gene annotations. |
| `biofilter_build` | *(required if annotating)* | Genome build to pass to Biofilter. Use `"19"` for GRCh37/hg19 or `"38"` for GRCh38/hg38. |
| `biofilter_close_dist` | `5E4` | Distance in base pairs defining "close" vs. "far" for nearest-gene annotation. A variant within this distance of a gene is labeled "close"; otherwise "far". Default is 50,000 bp (50 kb). |

---

## Infrastructure (Advanced)

These parameters set executable paths and environment-specific behavior. The defaults are correct for the SAIGE Docker/Singularity container and rarely need to change.

### Advanced

Only change these if using a custom container or non-standard compute environment.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `my_python` | `/opt/conda/bin/python` | Path to the Python executable. The default is correct for the SAIGE Docker/Singularity container. |
| `my_bgenix` | `/opt/conda/bin/bgenix` | Path to the `bgenix` executable. The default is correct for the SAIGE container. Only needed when `ftype = "BGEN"`. |
| `host` | `""` | Execution environment hint. Set to `"AOU"` for All of Us, `"DNAnexus"` for DNAnexus, `"LPC"` for the Penn LSF cluster, or leave empty for standard Docker/Singularity. Affects memory allocation and some process branching. |
| `gene_location_file` | `/app/NCBI.gene.loc` | Path to a gene coordinate file used for Manhattan plot labeling. The default is bundled in the container. |

---

# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"
    my_bgenix = "/opt/conda/bin/bgenix"

    //setting file type for step 2 (PLINK1/PLINK2/BGEN)
    //ftype = "PLINK1"
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
    is_imputed_data="TRUE" 
    minInfo=0.3 
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
## Advanced Nextflow Users: Take/Emit Info

### Input Channel (take) Description


NONE
### Output Channel (emit) Description


Singles_merge_output- a Channel of two paths to the merged raw and filtered output from SAIGE Step 2
Pheno_table- a Channel for a file that contains sample informations on each phenotype
