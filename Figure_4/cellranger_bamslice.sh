#!/bin/bash
#SBATCH --job-name bamslice
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 82g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/Software/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH


cellranger-dna bamslice --id=BAMs_3E5CFL \
                       --csv=3E5CFL_config.csv \
                       --bam=/home/regnerm/ENDO_OVAR_PROJ/ovar_3E5CFL_ATAC/possorted_bam.bam
cellranger-dna bamslice --id=BAMs_3BAE2L \
                       --csv=3BAE2L_config.csv \
                       --bam=/home/regnerm/ENDO_OVAR_PROJ/ovar_3BAE2L_ATAC/possorted_bam.bam
