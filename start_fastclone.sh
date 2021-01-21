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

python sequenza2pyclone.py \
/scratch/eknodel/Cancer_Genomics/01_somatic_mutation_calling/gatk_mutect2/2.somatic.filtered.pass.vcf.forpyclone \
/scratch/eknodel/FastClone_GuanLab/sequenza/2_segments.txt \
Patient2 \
/scratch/eknodel/FastClone_GuanLab/sequenza/

 
fastclone load-pyclone prop Patient2_sequenza2pyclone.txt  None solve ./Patient2_fastclone
