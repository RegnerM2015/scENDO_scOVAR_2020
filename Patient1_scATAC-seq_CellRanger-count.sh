#!/usr/bin/env bash

#SBATCH --job-name 3533EL-ATAC_A3
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3533EL-ATAC_A3.cellranger-count.job.update.out 
#SBATCH --error 3533EL-ATAC_A3.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=3533EL-ATAC_A3_Update \
                      --fastqs=./fastq_path2/H333JBGXB/3533EL-ATAC_A3 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=3533EL-ATAC_A3 \
                      --localcores=16
