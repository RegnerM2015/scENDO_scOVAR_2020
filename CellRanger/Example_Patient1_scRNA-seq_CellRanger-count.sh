#!/usr/bin/env bash

#SBATCH --job-name 3533EL-RNA_F6
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3533EL-RNA_F6.cellranger-count.job.update.out 
#SBATCH --error 3533EL-RNA_F6.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=3533EL-RNA_F6_Update \
                      --fastqs=./fastq_path/HMGVCBGX9/3533EL-RNA_F6 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=3533EL-RNA_F6 \
                      --localcores=16 \
                      --localmem=80
