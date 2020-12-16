#!/usr/bin/env bash

#SBATCH --job-name 3571DL-RNA_G1
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3571DL-RNA_G1.cellranger-count.job.update.out 
#SBATCH --error 3571DL-RNA_G1.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=3571DL-RNA_G1_Update \
                      --fastqs=./fastq_path/HMGVCBGX9/3571DL-RNA_G1 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=3571DL-RNA_G1 \
                      --localcores=16 \
                      --localmem=80
