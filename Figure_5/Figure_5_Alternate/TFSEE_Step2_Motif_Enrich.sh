#!/bin/bash
#SBATCH --job-name Motif_Select
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH

dir=/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/TFSEE_hg38_Feb_Select_Malignant


for i in *-updated.bed
do 

# Get sequence of marker enhancers for that cluster
bedtools getfasta -fi ${dir}/genome.fa -fo enhancers_${i}.fa -bed ${i}
	
meme enhancers_${i}.fa -dna -mod zoops -nmotifs 15 -minw 8 -maxw 15 -revcomp -oc enhancers_${i}_meme -maxsize 20000000

tomtom -evalue -thresh 10 -oc enhancers_${i}_tomtom ./enhancers_${i}_meme/meme.txt ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme

done 	



















