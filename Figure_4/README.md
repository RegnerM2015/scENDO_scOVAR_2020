
## Figure 4 plotting 

The starting inputs for this script are 1) the multi-sample (HGSOC cohort) Seurat object, 2) the multi-Sample ArchR project with HGSOC peak-to-gene associations, and 3) the HGSOC peak-to-gene link metadata. The outputs are the UMAP plots, proportion bar charts, expression boxplots for the gene of interest, and the ATAC browser track for the gene of interest presented in Figure 4. 

### [/Figure_4/Figure_4.R](https://github.com/RegnerM2015/scENDO_scOVAR_2020/tree/main/Figure_4)

1. Plot scRNA-seq/scATAC-seq UMAP plots colored by sample
1. Plot scRNA-seq/scATAC-seq UMAP plots colored by cell type
1. Plot proportion bar charts for scRNA-seq and scATAC-seq showing the contribution of each patient to each cell type subcluster
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in box plot


## Transcription factor motif analysis at distal-regulatory elements

We next wanted to predict TF occupancy at these putative regulatory elements for the gene of interest. To account for single-nucleotide variants in our motif analysis, we first created separate bam files for the maligant fractions of Patient 8 and 9. In other words, we wanted to only investigate the ATAC reads that came from cells in the predicted malignant clusters. We performed this operation using `bamslice` from cellranger-dna:

```
#!/bin/bash
#SBATCH --job-name bamslice
#SBATCH --cpus-per-task 16
#SBATCH -c 16
#SBATCH --mem 82g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/Software/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH


cellranger-dna bamslice --id=BAMs_3E5CFL \
                       --csv=3E5CFL_config.csv \
                       --bam=/home/regnerm/ENDO_OVAR_PROJ/ovar_3E5CFL_ATAC/possorted_bam.bam
cellranger-dna bamslice --id=BAMs_3BAE2L \
                       --csv=3BAE2L_config.csv \
                       --bam=/home/regnerm/ENDO_OVAR_PROJ/ovar_3BAE2L_ATAC/possorted_bam.bam
```
We used the malignant-specific bam file for Patient 9 (aka Patient ID 3E5CFL) as input to into variant calling followed by FIMO motif scanning in the genomic regions of interest:

```
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
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_enhancer_1.fa -bed enhancer_1.bed
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_enhancer_2.fa -bed enhancer_2.bed
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_enhancer_3.fa -bed enhancer_3.bed
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_enhancer_4.fa -bed enhancer_4.bed
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_enhancer_5.fa -bed enhancer_5.bed
bedtools getfasta -fi consensus.fa -fo GeneOfInterest_promoter.fa -bed promoter.bed

# Run FIMO motif scanning
fimo --bgfile motif-file --oc enhancer_1_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_enhancer_1.fa
fimo --bgfile motif-file --oc enhancer_2_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_enhancer_2.fa
fimo --bgfile motif-file --oc enhancer_3_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_enhancer_3.fa
fimo  --bgfile motif-file --oc enhancer_4_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_enhancer_4.fa
fimo --bgfile motif-file --oc enhancer_5_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_enhancer_5.fa
fimo --bgfile motif-file --oc promoter_fimo ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme GeneOfInterest_promoter.fa
```


Finally, we read the FIMO motif scanning results into R and screen for TF matches with a q-value < 0.10. To further rank these TF motif predictions, we rank them by descending order in normalized pseudo-bulk TF gene expression within the malignant fraction of Patient 9 and plot their corresponding expression patterns in box plots. Note that "0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","16-Epithelial cell" represent the malignant fraction of Patient 9. 

```
library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)


ovar_HGSOC_scRNA_processed <- readRDS("../ovar_HGSOC_scRNA_processed.rds")

labels <- c("0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","16-Epithelial cell")


hgsoc <- ovar_HGSOC_scRNA_processed[,ovar_HGSOC_scRNA_processed$cell.type %in% labels] 

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(as.matrix(counts)))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.files(pattern = "_fimo")

for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

# Make boxplot for enhancer 2 TF expression levels:
enh.2.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[12582,],
                          TF ="A.KLF6")
print(rownames(hgsoc)[12582])
enh.2.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                          TF="B.YY1")
print(rownames(hgsoc)[15548])
enh.2.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6857,],
                          TF="C.TFAP2A")
print(rownames(hgsoc)[6857])


# Make boxplot for enhancer 4 TF expression levels:
enh.4.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[9985,],
                          TF ="D.CEBPD")
print(rownames(hgsoc)[9985])
enh.4.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                          TF="E.YY1")
print(rownames(hgsoc)[15548])
enh.4.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6857,],
                          TF="F.TFAP2A")
print(rownames(hgsoc)[6857])

# Make boxplot for promoter TF expression levels:
prom.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[12582,],
                         TF ="G.KLF6")
print(rownames(hgsoc)[12582])
prom.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15260,],
                         TF="H.HIF1A")
print(rownames(hgsoc)[15260])
prom.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6925,],
                         TF="I.SOX4")
print(rownames(hgsoc)[6925])
prom.tfs.4 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                         TF="J.YY1")
print(rownames(hgsoc)[15548])


tfs <- rbind(enh.2.tfs.1,enh.2.tfs.2,enh.2.tfs.3,
                  enh.4.tfs.1,enh.4.tfs.2,enh.4.tfs.3,
                  prom.tfs.1,prom.tfs.2,prom.tfs.3,prom.tfs.4)
ggplot(tfs,aes(x=fct_rev(TF),y=Data))+
  geom_boxplot(lwd=0.45,outlier.size = 0.95,fatten = 0.95)+theme_classic()+ylim(c(0,4))+
  coord_flip()+
  ggsave("Vln.pdf",width =4,height = 8)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```
Interested in more exciting research in cancer genomics? Visit https://www.thefrancolab.org/ to learn more!
