#!/bin/bash
#SBATCH --job-name Motif 
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/Software/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH

dir=/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/TFSEE_Methods

for i in *_config-updated.csv
do
name=$(echo $i| cut -d'_' -f 2)
echo $name
echo $i
cellranger-dna bamslice --id=BAMs_${name} \
                       --csv=${i} \
                       --bam=/home/regnerm/ENDO_OVAR_PROJ/*${name}_ATAC/possorted_bam.bam
done 
for i in BAMs_*
do 
cd ${dir}/${i}/outs/subsets
for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done
	for j in *.bam
	do 
	# call variants
	bcftools mpileup -Ou -f ${dir}/genome.fa ${j} | bcftools call -mv -Oz -o variants_${j}.vcf.gz 
	
	bcftools index variants_${j}.vcf.gz

	# Apply cluster variants to output cluster reference genome
	cat ${dir}/genome.fa | bcftools consensus variants_${j}.vcf.gz > consensus_${j}.fa
	
	k=${j//.bam/}
	k="${k:7}"
	# Get sequence of marker enhancers (ranked list of enhancers sorted by FDR)
	bedtools getfasta -fi consensus_${j}.fa -fo enhancers_${k}.fa -bed ${dir}/Marker_Enhancers_ArchR_${k}-updated.bed
	
	meme enhancers_${k}.fa -dna -mod zoops -nmotifs 15 -minw 8 -maxw 15 -revcomp -oc enhancers_${k}_meme -maxsize 20000000


	tomtom -evalue -thresh 10 -oc enhancers_${k}_tomtom ./enhancers_${k}_meme/meme.txt ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme

	done
done 	



















