#!/bin/bash
#SBATCH --job-name Motif_Search 
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 32g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH

dir=/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/Figure_Making/Figure_4

# # call variants in the epithelial fraction of Patient 9
# The Patient_9 epithelial bam file was generated with Cell Rangerf's bamslice on the original Patient 9 bam file 
# (https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/using/bamslice)
bcftools mpileup -d 1000 -Ou -f /datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa ./BAMs_3E5CFL/outs/subsets/3E5CFL_Epithelial_cell.bam | bcftools call -mv -Oz -o variants.vcf.gz

bcftools index variants.vcf.gz
bcftools view --types snps variants.vcf.gz > variants_new.vcf
bgzip -c variants_new.vcf > variants_new.vcf.gz
tabix -p vcf variants_new.vcf.gz

#Apply cluster variants to output cluster reference genome
cat /datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa | bcftools consensus variants_new.vcf.gz > consensus.fa

# Get sequence of enhancers and promoter
bedtools getfasta -fi consensus.fa -fo LAPTM4B_enhancer_1.fa -bed enhancer_1.bed
bedtools getfasta -fi consensus.fa -fo LAPTM4B_enhancer_2.fa -bed enhancer_2.bed
bedtools getfasta -fi consensus.fa -fo LAPTM4B_enhancer_3.fa -bed enhancer_3.bed
bedtools getfasta -fi consensus.fa -fo LAPTM4B_enhancer_4.fa -bed enhancer_4.bed
bedtools getfasta -fi consensus.fa -fo LAPTM4B_enhancer_5.fa -bed enhancer_5.bed
bedtools getfasta -fi consensus.fa -fo LAPTM4B_promoter.fa -bed promoter.bed

# Run FIMO motif scanning
fimo --bgfile motif-file --oc enhancer_1_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_enhancer_1.fa
fimo --bgfile motif-file --oc enhancer_2_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_enhancer_2.fa
fimo --bgfile motif-file --oc enhancer_3_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_enhancer_3.fa
fimo  --bgfile motif-file --oc enhancer_4_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_enhancer_4.fa
fimo --bgfile motif-file --oc enhancer_5_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_enhancer_5.fa
fimo --bgfile motif-file --oc promoter_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme LAPTM4B_promoter.fa
