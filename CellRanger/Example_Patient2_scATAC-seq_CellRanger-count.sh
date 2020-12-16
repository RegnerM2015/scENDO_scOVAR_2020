#!/usr/bin/env bash

#SBATCH --job-name 3571DL-ATAC_A4
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3571DL-ATAC_A4.cellranger-count.job.update.out 
#SBATCH --error 3571DL-ATAC_A4.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=3571DL-ATAC_A4_Update \
                      --fastqs=./fastq_path2/H333JBGXB/3571DL-ATAC_A4 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=3571DL-ATAC_A4 \
                      --localcores=16
