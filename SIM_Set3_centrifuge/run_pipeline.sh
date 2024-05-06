#!/bin/sh
#SBATCH --job-name=centrifuge_SIM_Set3
#SBATCH --account=co_moorjani
#SBATCH --partition=savio4_htc 
#SBATCH --time=72:00:00 --mem 150gb 
#SBATCH --cpus-per-task=4

# Assuming you are in the directory containing the specific folders
module load samtools
module load r
module load gsl
module load bwa 
module load r/4.3.2
# r already downloaded in my user folder 
# pysam already downloaded in my user folder  

snakemake --core 1 -s ../rules/centrifuge.smk 