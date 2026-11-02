Documentation for SAIGE Pipelines
# Overview
These nextflow-enabled family of workflows are designed to allow users to quickly and easily perform association tests using [SAIGE](https://saigegit.github.io/SAIGE-doc/). These pipelines work with multiple types of genetic data, such as Common Variants, Rare Variants, and Gene Burden. Users provide input files and specify paramaters, after which the pipeline will automate preprocessing, step 1, step 2, and post processing & visualization.

Here are some SAIGE-specific resources:
- [Original SAIGE Paper for Your Reference](https://www.nature.com/articles/s41588-018-0184-y)
- [SAIGE-GENE Paper for Your Reference](https://www.nature.com/articles/s41588-022-01178-w)
- [SAIGE Documentation](https://saigegit.github.io/SAIGE-doc/)

## Software Requirements:
- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- Some text editor for updating the `nextflow.config` profiles such as `vim` or `nano`
- [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)
- [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
- [JDK version 11.0.5](https://www.oracle.com/java/technologies/javase/jdk11-archive-downloads.html)

## Overview of Steps to get started with one of our SAIGE Pipelines:
1. Clone the git repo and install dependencies. 
2. Set up a `nextflow.config` file with a profile for your compute system. An example for this can be found here: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/blob/main/Example_Configs/nextflow.config).
3. Create a pipeline-specific `.config` file specifying your parameters and input files. Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
4. Run the workflow using a command like this: `nextflow run /path/to/pmbb-nf-toolkit-{pipeline}/{pipeline}.nf`. More details can be found in module-specific READMEs 
     - [SAIGE GWAS README](https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family/blob/main/READMEs/SAIGE_GWAS_docs.md)
     - [SAIGE GATE README](https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family/blob/main/READMEs/SAIGE-GATE_docs.md)
     - [SAIGE ExWAS README](https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family/blob/main/READMEs/SAIGE_ExWAS_docs.md)
     - [SAIGE PheWAS - Gene Burden README](https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family/blob/main/READMEs/Single-Gene_Burden_PheWAS_docs.md)
     - [SAIGE PheWAS - Single Variant README](https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family/blob/main/READMEs/Single-Variant_PheWAS_docs.md)

# Detailed Steps for Runnning One of our Pipelines

Note: test data were obtained from the [SAIGE github repo](https://github.com/saigegit/SAIGE).

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
git clone https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family.git
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
```

3. Build the `saige.sif` singularity image
- you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
- this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly. 

```sh
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
singularity build saige.sif docker://pennbiobank/saige:latest
```

## Part II: Configure your run

1. Make a separate analysis/run/working directory.
   - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses. 
   - ❗This step is optional, but We Highly recommend making a  `tools` directory separate from your `run` directory. The only items that need to be in the run directory are the `nextflow.conf` file and the `${workflow}.conf` file.

```sh
WDIR="/path/to/analysis/run1"
mkdir -p 
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile ${MY_PROFILE}`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and symlinking it to your run directory.

3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
   - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters. 
   - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
   - you can compartamentalize your config file as much as you like by passing 
   - There are 2 ways to specify the config file during a run:
      - with the `-c` option on the command line: `nextflow run -c /path/to/workflow.conf`
      - in the `nextflow.conf`: at the top of the file add: `includeConfig '/path/to/workflow.conf'` 

## Part III: Run your analysis

- ❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly. 
- ❗We HIGHLY recommend doing a test run with the included test data in `${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/test_data`
- in the `test_data/` directory for each pipeline, we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -stub
# run an exwas for real
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf
# resume an exwas run if it was interrupted or ran into an error
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -resume
```

# How Compute Scales

The number of cohorts, phenotypes and chromosomes you choose sets how many jobs the pipeline submits. Other parameters don't add jobs but make each job run longer. This section explains both, so you can estimate the cost of a run before you start it. Each workflow README has a **Compute & Scaling** section with its own fan-out table and example.

## Where the compute goes

| Step | What it does | Run by the pipeline? | Runs once per | Relative cost |
|---|---|---|---|---|
| Step 0: sparse GRM (SAIGE `createSparseGRM`) | Builds a sparse relatedness matrix | **No.** Build it once beforehand and pass it with `step1_sparse_grm` (only used when `use_sparse_GRM = true`) | genotype dataset, reused by every cohort and phenotype | high, but one-time |
| Step 1 QC (`plink_qc_for_step1`) | Prunes the Step 1 genotypes and selects up to 150k GRM markers | Only when `use_sparse_GRM = false` | cohort | low to medium |
| Step 1 null model | Fits the null GLMM for one phenotype | Yes | cohort × phenotype | **highest** |
| Step 2 association tests | Tests variants or genes against the null model | Yes | cohort × phenotype × chromosome (or chunk) | medium to high, depending on the number of variants or genes |
| Merging, plots, report | Combines and summarizes results | Yes | varies by workflow | low |

### Dense vs sparse GRM

`use_sparse_GRM` doesn't change the number of jobs, but it changes the cost of every Step 1 job a lot.

- **Dense GRM (`use_sparse_GRM = false`)**: Step 1 QC picks the GRM markers once per cohort, but SAIGE recomputes the full GRM from those markers **inside every Step 1 job**. The GRM is never saved, so its cost is paid again for each phenotype. `max_vars_for_GRM` and the cohort size drive this cost.
- **Sparse GRM (`use_sparse_GRM = true`)**: Step 1 QC is skipped, and every Step 1 job reads the same pre-built sparse GRM. Each job is much cheaper. The one-time Step 0 cost is outside the pipeline.

## What adds jobs vs what slows each job

**Adds jobs.** These multiply the number of Step 1 and/or Step 2 jobs:
- `cohort_list`: every cohort gets its own Step 1 QC, Step 1 and Step 2 jobs. Sex-stratified cohorts are separate entries in `cohort_list`, so listing `EUR_ALL`, `EUR_M` and `EUR_F` triples the work.
- Phenotype lists (`bin_pheno_list`, `quant_pheno_list`, `survival_pheno_list`): one Step 1 job per cohort × phenotype. Only pairs that meet `min_bin_cases` / `min_quant_n` / `min_survival_cases` run, and sex-specific phenotypes only run in sex-stratified cohorts.
- `chromosome_list`: multiplies Step 2 only. Step 1 is never repeated per chromosome.
- Chunking (`enable_chunking`, GWAS and GATE): splits each chromosome's Step 2 job into chunks. That gives more, smaller jobs with about the same total CPU time, but more parallelism.
- Exposures (Single-Variant PheWAS exposure mode): one set of Step 1 jobs per distinct exposure sample set in each cohort.

**Slows each job.** These leave the job count the same:
- Samples per cohort: Step 1 and Step 2.
- `max_vars_for_GRM`: dense-GRM Step 1.
- `LOCO = "TRUE"` (GWAS, GATE, Single-variant PheWAS): Step 1 fits leave-one-chromosome-out predictions. ExWAS and Gene-burden PheWAS always fit Step 1 without LOCO, because SAIGE turns LOCO off when it estimates variance ratios by MAC category.
- `grouptest_maf` and `grouptest_annotation` (ExWAS, Gene-burden PheWAS): every MAF cutoff × annotation mask is tested inside the same Step 2 job, so 3 cutoffs × 5 masks means longer jobs, not more jobs.
- `burden_only = false`: SKAT and SKAT-O tests are slower than burden-only tests.
- `use_firth`: extra fitting in binary Step 2 jobs, only for tests with p < `firth_cutoff`.
- Variants or genes per chromosome: Step 2. `snplist` and `gene_list_file` in the PheWAS workflows shrink this.

## Job counts by workflow

C = cohorts, P = cohort × phenotype pairs that pass the filters above (per cohort), K = chromosomes in `chromosome_list` (or chunks with `enable_chunking`), S = distinct exposure sample sets per cohort.

| Workflow | Run with | Step 1 jobs | Step 2 jobs |
|---|---|---|---|
| [SAIGE GWAS](READMEs/SAIGE_GWAS_docs.md#compute--scaling) | `saige_gwas.nf` | C × P | C × P × K |
| [SAIGE GATE](READMEs/SAIGE-GATE_docs.md#compute--scaling) | `saige_gwas.nf` with `survival_pheno_list` | C × P | C × P × K |
| [SAIGE ExWAS](READMEs/SAIGE_ExWAS_docs.md#compute--scaling) | `saige_exwas.nf` | C × P | C × P × K |
| [Gene-burden PheWAS](READMEs/Single-Gene_Burden_PheWAS_docs.md#compute--scaling) | `saige_gene_phewas.nf` | C × P | C × P × K |
| [Single-variant PheWAS](READMEs/Single-Variant_PheWAS_docs.md#compute--scaling) | `saige_variant_phewas.nf` | C × P | C × P × K |
| [Phenotype exposure PheWAS](READMEs/Single-Variant_PheWAS_docs.md#compute--scaling) | `saige_variant_phewas.nf` with exposures | C × S × P | C × S × P |

Step 1 counts are the same for every genetic-exposure workflow: testing 5 genes or 5 SNPs costs the same Step 1 as a genome-wide scan. The difference between workflows is in Step 2, i.e. how much each job tests.

**Example:** 2 cohorts × 100 phenotypes × 22 chromosomes is 200 Step 1 jobs and 4,400 Step 2 jobs. Adding a third cohort adds 100 Step 1 and 2,200 Step 2 jobs. Adding a second MAF cutoff adds none.