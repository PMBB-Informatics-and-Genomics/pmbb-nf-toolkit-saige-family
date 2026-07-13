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