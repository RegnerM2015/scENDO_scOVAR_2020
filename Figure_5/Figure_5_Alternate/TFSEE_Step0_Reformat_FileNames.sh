#!/bin/bash
#SBATCH --job-name Motif 
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes

# Remove whitespace from file names 

#for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done


# Concatenate all marker enhancer BED files into one 

cat Marker_Enhancers_ArchR* > Marker_Enhancers_ArchR.bed
sort Marker_Enhancers_ArchR.bed | uniq -u > Marker_Enhancers_ArchR-nodups.bed


