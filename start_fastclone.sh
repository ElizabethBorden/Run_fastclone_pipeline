#!/bin/bash
#SBATCH --job-name=QC  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 168:00:00
#SBATCH --mem=40000

newgrp combinedlab

source activate var_call_env

#export
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile fastclone_pipeline-snakemake.py -j 15 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 8 -t 1:00:00"
