#!/bin/bash
#SBATCH --job-name EEC_ATAC
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 92g
#SBATCH --partition allnodes

#export R_LIBS_USER="/home/regnerm/R/x86_64-pc-linux-gnu-library/4.0:${R_LIBS_USER}"
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate scENDO_scOVAR 
R CMD BATCH './EmpFDR_Estimation-Patient2.R'







