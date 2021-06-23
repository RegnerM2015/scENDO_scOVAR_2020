#!/bin/bash
#SBATCH --job-name Reformat 
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 2g
#SBATCH --partition allnodes

# Remove whitespace from file names 

for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done
