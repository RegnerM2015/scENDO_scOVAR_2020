#!/usr/bin/env bash

#SBATCH --job-name HMGVCBGX9
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HMGVCBGX9_demultiplex.job.out
#SBATCH --error HMGVCBGX9_demultiplex.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scRNA-seq_Endometrial.05.15.2019

cellranger mkfastq --id=HMGVCBGX9 \
                   --run=${DATA}/190514_NS500270_0297_AHMGVCBGX9 \
                   --csv=${OUT}/Samplesheet.csv \
                   --qc \
                   --localcores=16
