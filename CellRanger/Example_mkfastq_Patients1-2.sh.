#!/usr/bin/env bash

#SBATCH --job-name H333JBGXB_2
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output H333JBGXB_demultiplex.2.job.out
#SBATCH --error H333JBGXB_demultiplex.2.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scATAC-seq_Endometrial.05.16.2019

cellranger-atac mkfastq --id=H333JBGXB_2 \
                        --run=${DATA}/190515_NS500270_0298_AH333JBGXB \
                        --csv=${OUT}/Samplesheet.2.csv \
                        --qc \
                        --localcores=16
