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
