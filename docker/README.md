# Docker Container / Singularity Image README

Purpose of SAIGE Docker Container / Singularity Image: provides all the correct languages, packages, and scripts needed for this pipeline in a way that nextflow can easily access.

## Relevant Container Locations
* Docker container on Docker Hub: <code>pennbiobank/saige:latest</code>
* Docker container on Google Container Registry: <code>gcr.io/ritchie-aou-psom-9015/saige:latest</code>
* Singularity image: <code>saige.sif</code>

## Using the SAIGE Docker Container / Singularity Image in the Nextflow Pipeline
* The names and locations of these images are written in the <code>nextflow.config</code> file.
* The pipeline will either use a Docker container or Singularity image depending on which software the platform the pipeline is being run on supports.
* LPC (cluster) supports Singularity, and DNA Nexus and All of Us support Docker.
* There are three profiles written in the <code>nextflow.config</code> file that can be used to run the pipeline on LPC (cluster), DNA Nexus, or All of Us. These profiles specify paths to the Docker container or Singularity image each platform supports.
* When running the pipeline on different platforms, change the <code>-profile</code> flag in the nextflow command. <code>nextflow run saige_exwas.nf -resume -profile cluster</code> will work on LPC (cluster), <code>nextflow run saige_exwas.nf -resume -profile standard</code> will work on DNA Nexus, and <code>nextflow run saige_exwas.nf -resume -profile all_of_us</code> will work on All of Us.
* TBD: When using a Docker container, the pipeline pulls Docker container from either Dockerhub or Google Container Registry, depending on the platform this pipeline is being run on (thus, you do not need to pull the Docker image for this pipeline to work)
* The singularity image must be in the directory the nextflow command is being called from.
* To use the Singularity image on LPC, you must do <code>module load singularity</code> first (loads the Singularity software).

## Commands Used to Build Docker Image
* This Docker image was built on a Mac with an M1 chip, so the command was modified for that machine (but can be used on any machine).
* Dockerfile must be in the current directory the build command was called from.
### Building Docker Image and Pushing to my Dockerhub Account
<code>docker buildx build --platform  linux/amd64 --push -t pennbiobank/saige:latest .</code>
### Building Docker Image and Pushing to Public Google Container Registry
<code>docker buildx build --platform  linux/amd64 --push -t gcr.io/ritchie-aou-psom-9015/saige:latest .</code>

## Command Used to Build Singularity Image
* Useful for systems that do not support Docker.
* This command pulls Docker container from a publicly available site (Dockerhub or Google Container Registry) and converts it to a Singularity image locally.
<code>singularity build saige.sif docker://pennbiobank/saige:latest</code>

## Commands Used to Pull Docker Image
* From Dockerhub: <code>docker pull pennbiobank/saige:latest</code>
* From Google Container Registry: <code>docker pull gcr.io/ritchie-aou-psom-9015/saige:latest</code>
