#!/bin/bash
#SBATCH --job-name Intersect
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 64g
#SBATCH --partition allnodes

#export R_LIBS_USER="/home/regnerm/R/x86_64-pc-linux-gnu-library/4.0:${R_LIBS_USER}"
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate scENDO_scOVAR 
R CMD BATCH './Intersect_DiffPeaks_P2Gs.R'







