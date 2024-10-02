#!/bin/bash

#Chris Carson
#April 12, 2024
# updated: Zach Rodriguez June 2024

# I got tired of doing all the symlinks by hand so I made this script to automate it.
# Just give it the name of the folder you want to make
# when running this "bash makesaigefolder.sh *insertfoldernamehere* */full/path/to/geno_pheno_workbench/*"

target=$1
gpw=$2

mkdir -p $target
cd $target

# symlink to keep original
ln -s ${gpw}/SAIGE_FAMILY/processes ${target}/processes
ln -s ${gpw}/SAIGE_FAMILY/scripts ${target}/scripts
ln -s ${gpw}/SAIGE_FAMILY/workflows ${target}/workflows
ln -s ${gpw}/SAIGE_FAMILY/READMES ${target}/READMES
ln -s ${gpw}/SAIGE_FAMILY/saige.sif ${target}/saige.sif
ln -s ${gpw}/SAIGE_FAMILY_GPU/saige-doe-3.sif ${target}/saige-doe-3.sif

# copy these for editing
cp -r ${gpw}/SAIGE_FAMILY/configs ${target}/configs
cp -r ${gpw}/SAIGE_FAMILY/nextflow.config ${target}/nextflow.config
cp -r ${gpw}/nextflow_gpu.config ${target}/nextflow_gpu.config
