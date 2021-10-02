We adapted our method called Total Functional Score of Enhancer Elements (TFSEE) to predict which TFs are enriched at active enhancer-like elements within malignant cell types (Franco et al., 2018; Malladi et al., 2020). To do this, we leveraged our matched scRNA-seq and scATAC-seq data to infer 1) TF expression, 2) enhancer-like element activity and 3) TF motif prediction within enhancer-like elements for each malignant population of interest.  TFSEE allows for the simultaneous assessment of enhancer activity and TF gene expression, and identification of TFs that drive gene expression programs in malignant cell states. Referring back to the entire patient cohort, 11 malignant cell type subclusters were chosen for TFSEE analysis based on patient specificity, inferred CNV events and/or cancer biomarker expression patterns.
<p align="center">
Adapted from Figure 6A
</p>
<p align="center">
<img src="https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/TFSEE_cartoon-update.PNG" width=775" height="250">
</p>


## Part I: set intersection of distal regulatory elements with differentially accessible peaks 

We first identify which of the distal regulatory elements (correlation >= 0.45 and p-value <= 1e-12) are enriched with statistical significance in the malignant cell type populations from the full cohort. The set intersection, or character intersection, of the peak names are written out to bed files for each malignant cell type included in the TFSEE analysis. After finding the distal regulatory elements that are enriched with statistical significance in each malignant cluster, we concatenate all of the results and remove duplicate peak names (since one peak could be enriched in multiple clusters). This concatenated list of differentially accessible distal regulatory elements is written out to a bed file. These operations are performed in the following script which is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/Intersect_DiffPeaks_P2Gs.R

```
library(ArchR)
library(ChIPpeakAnno)
library(stats)
ArchR::addArchRThreads(threads = 32)
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./getMatrixFromProject.R")
h5disableFileLocking()
SAMPLE.ID <- "All"
key.word.1 <- "pithelia"
key.word.2 <- "-Ciliated"
proj <- readRDS("./final_archr_proj_archrGS.rds")

p2g.df.obs <- readRDS("./All_P2G_Observed.rds")
#Subset to postive correlation P2Gs
p2g.df.sub <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45& peakType == "Distal")


# Marker Peaks for malignant clusters with 100% patient specificity 
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList,"markerspeaks.rds")


print(names(markerList))

idx.1 <- grep(key.word.1,names(markerList))
idx.2 <- grep(key.word.2,names(markerList))
idx.3 <- grep("16-Fibroblast",names(markerList))
idx.4 <- grep("0-Fibroblast",names(markerList))
idx.5 <- grep("27-Fibroblast",names(markerList))

idx <- c(idx.1,idx.2,idx.3,idx.4,idx.5)

markerList.sub <- markerList[idx]

print(names(markerList.sub))

# Only tumor cell clusters that are 100% patient specific 

df <- as.data.frame(proj@cellColData) %>% dplyr::group_by(predictedGroup_ArchR) %>% 
  dplyr::count(Sample) 

idents <- table(df$predictedGroup_ArchR)
idents.sub <- idents[idents ==1 & names(idents) %in% names(markerList.sub)]

markerList.sub <- markerList.sub[names(idents.sub)] 
print(names(markerList.sub))
# Intersect marker peaks for each cluster with P2Gs and write to bed file
for ( i in names(markerList.sub)){
  data <- unlist(markerList.sub[i])
  markers <- paste0(data$seqnames,":",data$start,"-",data$end)
  markers.int <- intersect(p2g.df.sub$peakName,markers)
  markers.int <- stringr::str_replace(string=markers.int,pattern = ":",replacement = "-")
  bed <- data.frame(do.call(rbind, strsplit(markers.int, "-", fixed=TRUE)))
  write.table(bed,paste0("Marker_Enhancers_ArchR_",i,".bed"),row.names = F,col.names = F,quote = F,sep = "\t")
}

# Concatenate marker peaks from each cluster and remove redundant peaks
data <- unlist(markerList.sub)
markers <- paste0(data$seqnames,":",data$start,"-",data$end)
markers.int <- intersect(p2g.df.sub$peakName,markers)
markers.int <- unique(markers.int)
markers.int <- stringr::str_replace(string=markers.int,pattern = ":",replacement = "-")
bed <- data.frame(do.call(rbind, strsplit(markers.int, "-", fixed=TRUE)))
write.table(bed,"Marker_Enhancers_ArchR-nodups.bed",row.names = F,col.names = F,quote = F,sep = "\t")



writeLines(capture.output(sessionInfo()), "sessionInfo_Intersect_DiffPeaks_P2Gs.txt")
```

## Part II: Generating the input matrices for the TFSEE computation

For technical reasons associated with the syntax of cell type cluster labels, we must first remove spaces from all file names using the following script which is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/TFSEE_Step0_Reformat_FileNames.sh
```
#!/bin/bash
#SBATCH --job-name Reformat 
#SBATCH --cpus-per-task 1
#SBATCH -c 1
#SBATCH --mem 2g
#SBATCH --partition allnodes

# Remove whitespace from file names 

for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done
```

Next, we identify which of the malignant-enriched distal regulatory elements have non-zero counts across all 11 malignant cell type clusters. We do this because we do not have replicates of each malignant cell type cluster to distinguish between biological and technical zeros. We update all of the bed files to include only the distal regulatory elements that have non-zero counts across all 11 malignant cell type clusters. These operations are performed in the following script which is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/TFSEE_Step1_NonZero_Enhancers.R

```
library(ArchR)
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./getMatrixFromProject.R")
library(stringr)
library(utils)
library(dplyr)
ArchR::addArchRThreads(threads = 32)
h5disableFileLocking()
# ONlY USE 100% patient-specfic 
labels <- list.files(pattern = "Marker_Enhancers_ArchR")
labels <- str_remove(labels,"Marker_Enhancers_ArchR_")
labels <- str_remove(labels,".bed")
labels <- labels[-12]# Remove extra
print(labels)

enhancers <- read.delim("Marker_Enhancers_ArchR-nodups.bed",header = F)
enhancers <- paste0(enhancers$V1,":",enhancers$V2,"-",enhancers$V3)

# Pseudobulk ATAC enhancer matrix
proj.archr <- readRDS("./final_archr_proj_archrGS-P2Gs.rds")
peak.mat <- getMatrixFromProject.mod(proj.archr,useMatrix = "PeakMatrix")
cell.names <- colnames(assay(peak.mat))
peak.names <- paste0(seqnames(rowRanges(peak.mat)),":", start(ranges(rowRanges(peak.mat))),"-", end(ranges(rowRanges(peak.mat))))


peak.mat <- assay(peak.mat)

dim(peak.mat)
length(cell.names)
length(peak.names)

colnames(peak.mat) <- cell.names
rownames(peak.mat) <- peak.names
head(peak.mat[,1:4])


peaks.pseudobulk <- data.frame(rownames(peak.mat))
# Run twice to hit both whitespaces 
proj.archr$predictedGroup_ArchR <- str_replace(proj.archr$predictedGroup_ArchR," ","_")
proj.archr$predictedGroup_ArchR <- str_replace(proj.archr$predictedGroup_ArchR," ","_")
for (i in labels){
  cells <- rownames(dplyr::filter(as.data.frame(proj.archr@cellColData),predictedGroup_ArchR == i))
  
  peak.mat.sub <- peak.mat[,colnames(peak.mat) %in% cells]
  
  peaks.bulk <- Matrix::rowSums(peak.mat.sub)
  
  peaks.pseudobulk$i <- peaks.bulk
  colnames(peaks.pseudobulk)[dim(peaks.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}

rownames(peaks.pseudobulk) <- peaks.pseudobulk[,1]
peaks.pseudobulk <- peaks.pseudobulk[,-1]

print(labels)
print(colnames(peaks.pseudobulk))
dim(peaks.pseudobulk)
head(peaks.pseudobulk[,1:6])



peaks.pseudobulk <- peaks.pseudobulk[rownames(peaks.pseudobulk) %in% unique(enhancers),]
dim(peaks.pseudobulk)
# Remove enhancer peaks that have ANY zero counts in ANY sample (no replicates)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]
dim(peaks.pseudobulk)
head(peaks.pseudobulk[,1:6])


# Subset original Marker enhancer lists to new updated nonzero enhancers
for ( i in 1:length(labels)){
  file <- paste0("Marker_Enhancers_ArchR_",labels[i],".bed")
  
  bed <- read.delim(file,header = F)
  rownames(bed) <- paste0(bed$V1,":",bed$V2,"-",bed$V3)
  
  bed.new <- bed[rownames(bed) %in% rownames(peaks.pseudobulk),]
  
  write.table(bed.new,paste0("Marker_Enhancers_ArchR_",labels[i],"-updated.bed"),row.names = F,col.names = F,quote = F,sep = "\t")
}



writeLines(capture.output(sessionInfo()), "sessionInfo_TFSEE_Step1_NonZero_Enhancers.txt")
```

Using the bed files recording genomic coordinates of the distal regulatory elements enriched in each malignant cell type cluster, we conduct MEME motif discovery followed by TOMTOM motif matching. This step will eventually result in the TF motif prediction matrix required for the TFSEE computation. The following script performs these operations and is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/TFSEE_Step2_Motif_Enrich.sh

```
#!/bin/bash
#SBATCH --job-name Motif_All_Malign 
#SBATCH --cpus-per-task 32
#SBATCH -c 32
#SBATCH --mem 64g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH

dir=/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/Figure_Making/Figure_6


for i in *-updated.bed
do 
# FASTA: https://ondemand.bioinf.unc.edu/pun/sys/files/fs/datastore/lbcfs/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta
# Get sequence of marker enhancers for that cluster
bedtools getfasta -fi ${dir}/genome_hg38.fa -fo enhancers_${i}.fa -bed ${i}
	
meme enhancers_${i}.fa -dna -mod zoops -nmotifs 15 -minw 8 -maxw 15 -revcomp -oc enhancers_${i}_meme -maxsize 20000000

tomtom -evalue -thresh 10 -oc enhancers_${i}_tomtom ./enhancers_${i}_meme/meme.txt ${dir}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme

done 	
```


To process the MEME/TOMTOM outputs for each malignant cell type cluster, we use the following script adapted from https://git.biohpc.swmed.edu/gcrb/tfsee/-/blob/master/analysis/GRO_seq_TFSEE/matrix_analysis.py. Our adapted script results in the motif prediction matrix and is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/TFSEE_Step3_Motif_Pred_Matrix.py

```
import Bio.motifs
import pandas as pd
import numpy as np
import csv
import math
import re
import string
from sklearn import preprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.stats import pearsonr, spearmanr, cumfreq




### Parse MEME and TOMTOM Motif data

# Loop through meme output
meme_cell_dict = {"./enhancers_Marker_Enhancers_ArchR_0-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_0-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_10-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_10-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_16-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_16-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_17-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_17-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_19-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_19-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_21-Unciliated_epithelia_1-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_21-Unciliated_epithelia_1-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_27-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_27-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_3-Epithelial_cell-updated.bed_meme":"enhancers_Marker_Enhancers_ArchR_3-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_31-Unciliated_epithelia_1-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_31-Unciliated_epithelia_1-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_34-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_34-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_9-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_9-Epithelial_cell-updated.bed_tomtom"}


# Read Target ID to Motif into dictionary
motif_id_dict = {}
with open("./JASPAR_2020_Homo_Sapiens.txt", "r") as data:
    motif_ids = csv.DictReader(data, delimiter="\t")
    for line in motif_ids:
        motif_id_dict[line['ID']] = line['NAME']

meme_tomtom = pd.DataFrame()
for meme,tom in meme_cell_dict.items():
    # load meme output
    meme_file = '%s/meme.xml' % (meme)
    record = Bio.motifs.parse(open(meme_file),'meme')
    # Loop through all motifs and make dataframe
    meme_positions = pd.DataFrame()
    for motif in record:
        name = motif.name
        ones = [1] * len(motif.instances)
        names = []
        for instance in motif.instances:
            names.append(instance.sequence_name)
        new = pd.DataFrame({name: ones},index = names)
        temp = pd.concat([meme_positions, new], axis=1).fillna(0)
        meme_positions = temp
    # Read tomtom file
    tomtom_file = "%s/tomtom.txt" % (tom)
    tomtom_dict = {}
    with open(tomtom_file, "r") as data:
        tomtom = csv.DictReader(data, delimiter="\t")
        for line in tomtom:
            target = line['Target ID']
            motif = line['#Query ID']
            pval = float(line['p-value'])
            tfs = motif_id_dict[target].upper()
            motif_pvalue = { motif:  [pval]}
            # JASPAR :: means that any TF can either protein, split the protein
            tf_list = tfs.split("::")
            for tf in tf_list:
                # Reduce split form splice to single value [ID]_#
                single_isoform = tf.split("_")[0]
                if single_isoform in tomtom_dict.keys():
                    if motif in tomtom_dict[single_isoform].keys():
                        tomtom_dict[single_isoform][motif].append(pval)
                    else:
                        tomtom_dict[single_isoform].update(motif_pvalue)
                else:
                    tomtom_dict[single_isoform] = motif_pvalue
    # Make dataframe
    tomtom_motif = pd.DataFrame()
    for key,motif in tomtom_dict.items():
        pvalue_dict = {}
        # Loop through motifs to see if length greater than 1, if so do pvalue scaling
        for m,p in motif.items():
            if len(p) > 1:
                stouffer_statistic, stouffer_pval = scipy.stats.combine_pvalues(p,method = 'stouffer', weights = None)
                pvalue_dict[m] = stouffer_pval
            else:
                pvalue_dict[m] = p[0]
        pvalues = np.array(list(pvalue_dict.values()))
        new = pd.DataFrame({key: pvalues},index = pvalue_dict.keys())
        temp = pd.concat([tomtom_motif, new], axis=1).fillna(0).sort_index(level=int)
        tomtom_motif = temp
    # Reorder
    tomtom_motif_reorder = tomtom_motif.reindex( list(meme_positions.columns.values)).fillna(0)
    # dot product
    meme_tomtom_cell = meme_positions.dot(tomtom_motif_reorder)
    # Scale and add
    vals = meme_tomtom_cell.max(axis=0)
    print(vals)
    meme_tomtom_cell = np.log10(meme_tomtom_cell+0.01)
    meme_tomtom_cell = np.negative(meme_tomtom_cell)
    scaler = preprocessing.MinMaxScaler()
    meme_tomtom_cell_transform = meme_tomtom_cell.T
    print(meme_tomtom_cell_transform.shape)
    norm = scaler.fit_transform(meme_tomtom_cell_transform.values) # norm across enhancers for each enhancer
    meme_tomtom_cell_std = pd.DataFrame(data=norm.T, columns=list(meme_tomtom_cell.columns.values), index = meme_tomtom_cell.index )
    # Add to previous data
    temp = meme_tomtom.add(meme_tomtom_cell_std, fill_value=0).fillna(0).sort_index(level=int)
    meme_tomtom = temp



# Transform meme tom_tom
motif_enhancers = meme_tomtom.T

# Rename column headers
#motif_enhancers.rename(columns=lambda x: x.split('-')[0], inplace=True)
#motif_enhancers.rename(columns=lambda x: x.replace(':', "_"), inplace=True)

# Standardize to range 0-1
scaler = preprocessing.MinMaxScaler()
motif_enhancers_transform = motif_enhancers.T
print(motif_enhancers_transform.shape)
norm = scaler.fit_transform(motif_enhancers_transform.values) # norm across enhancers for each enhancer
motif_enhancers_scaled = pd.DataFrame(data=norm.T, columns=list(motif_enhancers.columns.values), index = motif_enhancers.index)



print(motif_enhancers.index)
motif_enhancers.to_csv("motif_enhancers.txt",sep = "\t")
motif_enhancers_scaled.to_csv("motif_enhancers_scaled.txt",sep = "\t")
```

## Part III: Assemble inputs and run TFSEE:

Finally, we assemble the TF motif prediction matrix, the TF expression matrix, and the enhancer matrix. To normalize the pseudo-bulk TF expression matrix and the pseudo-bulk peak/enhancer matrix, we apply the regularized logarithm transformation from DESeq2. Then we subset each matrix (so that they all have shared features) and scale from 0-1. The scaled enhancer activity matrix was multiplied with the scaled TF motif prediction matrix to form an intermediate matrix product. This matrix product was element-wise multiplied with the scaled TF expression matrix to form the final TFSEE matrix used in downstream analysis. These operations were performed in the following script which is also listed below: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/TFSEE_Step4_Matrix_Calculations.R


```
library(ArchR)
library(ComplexHeatmap)
library(DESeq2)
library(Seurat)
library(ggplot2)
library(dplyr)

###############################################
# Part 1: Read in TF motif prediction matrix 
###############################################
# Read in TF motif prediction matrix
enhancer.motifs<- read.delim("./motif_enhancers_scaled.txt",sep = "\t")
rownames(enhancer.motifs) <- enhancer.motifs[,1]
enhancer.motifs <- enhancer.motifs[,-1]
# Rename syntax of peaks
colnames(enhancer.motifs) <- sub("\\.",":",colnames(enhancer.motifs))
colnames(enhancer.motifs) <- sub("\\.","-",colnames(enhancer.motifs))


# Set labels of cell types of interest

# Only use clusters with 100% patient specificity 
labels <- c("31-Unciliated epithelia 1",
            "21-Unciliated epithelia 1",
            "19-Epithelial cell",
            "34-Epithelial cell",
            "3-Epithelial cell",
            "10-Epithelial cell",
            "16-Fibroblast",
            "17-Epithelial cell",
            "0-Fibroblast",
            "27-Fibroblast",
            "9-Epithelial cell")# Screen for "pure" patient clusters


###############################################
# Part 2: Pseudobulk TF expression 
###############################################

# Pseudobulk RNA TF expression
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")
rna.counts <- rna@assays$RNA@counts

rna.pseudobulk <- data.frame(rownames(rna))
for (i in labels){
  cells <- rownames(dplyr::filter(rna@meta.data,cell.type == i))
  
  rna.counts.sub <- rna.counts[,colnames(rna.counts) %in% cells]
  
  rna.counts.bulk <- rowSums(as.matrix(rna.counts.sub))
  
  rna.pseudobulk$i <- rna.counts.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
}
rownames(rna.pseudobulk) <- rna.pseudobulk[,1]
rna.pseudobulk <- rna.pseudobulk[,-1]
dim(rna.pseudobulk)

# Remove any features/rows that have zero counts in ANY sample
dim(rna.pseudobulk)
rna.pseudobulk[rna.pseudobulk == 0] <- NA
rna.pseudobulk <- rna.pseudobulk[complete.cases(rna.pseudobulk),]
dim(rna.pseudobulk)


# rlog normalize data
library(DESeq2)
dds.rna <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds.rna <- estimateSizeFactors(dds.rna)
sizeFactors(dds.rna)

rlog.rna <- rlog(dds.rna,blind = T)
saveRDS(dds.rna,"dds_rna.rds")
saveRDS(rlog.rna,"rlog_rna.rds")


###############################################
# Part 3: Pseudobulk enhancer activity
###############################################

# Pseudobulk ATAC enhancer matrix

atac <- readRDS("./final_archr_proj_archrGS.rds")

source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./getMatrixFromProject.R")
peaks.mat <- getMatrixFromProject.mod(atac,useMatrix = "PeakMatrix")
cell.names <- colnames(assay(peaks.mat))
peak.names <- paste0(seqnames(rowRanges(peaks.mat)),":", start(ranges(rowRanges(peaks.mat))),"-", end(ranges(rowRanges(peaks.mat))))


peaks.mat <- assay(peaks.mat)

dim(peaks.mat)
length(cell.names)
length(peak.names)

colnames(peaks.mat) <- cell.names
rownames(peaks.mat) <- peak.names
head(peaks.mat[,1:4])


# Construct pseudobulk peak/enhancer matrix
peaks.pseudobulk <- data.frame(rownames(peaks.mat))
for (i in labels){
  cells <- rownames(dplyr::filter(as.data.frame(atac@cellColData),predictedGroup_ArchR == i))
  
  peaks.sub <- peaks.mat[,colnames(peaks.mat) %in% cells]
  peaks.sub <- as(peaks.sub, "dgCMatrix")
  peaks.bulk <- Matrix::rowSums(peaks.sub)
  
  peaks.pseudobulk$i <- peaks.bulk
  colnames(peaks.pseudobulk)[dim(peaks.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}
rownames(peaks.pseudobulk) <- peaks.pseudobulk$rownames.peaks.
peaks.pseudobulk <- peaks.pseudobulk[,-1]

# Remove features/peaks that have zero counts in ANY sample
dim(peaks.pseudobulk)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]
dim(peaks.pseudobulk)

#DESeq ATAC

library(DESeq2)
dds.atac <- DESeqDataSetFromMatrix(countData = peaks.pseudobulk,
                              colData = data.frame(meta=colnames(peaks.pseudobulk)),design = ~ 1)

dds.atac <- estimateSizeFactors(dds.atac)
sizeFactors(dds.atac)

rlog.atac <- rlog(dds.atac,blind = T)

saveRDS(rlog.atac,"rlog_atac.rds")
saveRDS(dds.atac,"dds_atac.rds")
rlog.rna <- readRDS("rlog_rna.rds")
rlog.atac <- readRDS("rlog_atac.rds")

###############################################
# Part 4: Subsetting and scaling from 0 to 1
###############################################

# Subset and reorder motif matrix, peak matrix and expression matrix
enhancer.mat <- assay(rlog.atac)
enhancer.mat <- as.data.frame(enhancer.mat)
enhancer.mat <- enhancer.mat[rownames(enhancer.mat) %in% colnames(enhancer.motifs),]

# Check order
motif.mat <- enhancer.motifs
motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(enhancer.mat)]
motif.mat <- motif.mat[ ,order(match(colnames(motif.mat), rownames(enhancer.mat))) ]
length(which(colnames(motif.mat) == rownames(enhancer.mat)))

expr.mat <- assay(rlog.rna)
expr.mat <- as.data.frame(expr.mat)

# Check order
tfs <- intersect(rownames(motif.mat),rownames(expr.mat))
expr.mat <- expr.mat[rownames(expr.mat) %in% tfs,]
motif.mat <- motif.mat[rownames(motif.mat) %in% tfs,]
expr.mat <- expr.mat[order(match(rownames(expr.mat), rownames(motif.mat))), ]

length(which(rownames(expr.mat)==rownames(motif.mat)))

# Rescale peak matrix and rna expr matrix
library(scales)
for (i in 1:ncol(expr.mat)){
  expr.mat[,i] <- rescale(expr.mat[,i],to = c(0,1))
}
library(scales)# Scale cell types
for (i in 1:ncol(enhancer.mat)){
  enhancer.mat[,i] <- rescale(enhancer.mat[,i],to = c(0,1))
}


dim(enhancer.mat)
dim(motif.mat)
dim(expr.mat)



length(which(rownames(motif.mat) == rownames(expr.mat)))

length(which(rownames(enhancer.mat) == colnames(motif.mat)))



length(which(rownames(motif.mat) == rownames(expr.mat)))

length(which(rownames(enhancer.mat) == colnames(motif.mat)))

dim(enhancer.mat)
dim(motif.mat)
dim(expr.mat)


###############################################
# Part 5: TFSEE matrix operations
###############################################

# Peform matrix multiplication of enhancer activity by motif prediction
int <- t(enhancer.mat) %*% t(motif.mat)

dim(int)
dim(expr.mat)

length(which(colnames(int)==rownames(expr.mat)))



# Multiply intermediate matrix by TF expression matrix 
tfsee <- int*t(expr.mat)
dim(tfsee)
head(tfsee)

saveRDS(tfsee,"tfsee_matrix.rds")
saveRDS(int,"int_matrix.rds")
saveRDS(expr.mat,"expr_matrix.rds")
saveRDS(motif.mat,"motif_matrix.rds")

tfsee <- readRDS("tfsee_matrix.rds")
# Write table of TFs for canSAR search:
write.csv(data.frame(colnames(tfsee)),"TFs_for_canSAR.csv")

##################################################################
# END OF TFSEE COMPUTATION
##################################################################
```
## Part IV: Plotting Figure 6

We use the following script to plot the rank order frequency distribution plots and the heatmap shown in Figure 6: https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Figure_6/Figure_6_Plotting.R. Note that the canSAR druggability predictions were obtained by inputing the list of TFs into the CPAT tool at https://cansarblack.icr.ac.uk/cpat. 

```
library(ArchR)
library(ComplexHeatmap)
library(DESeq2)
library(Seurat)
library(ggplot2)
library(dplyr)

# Read in tfsee matrix and write TFs out to csv for canSAR
tfsee <- readRDS("./tfsee_matrix.rds")
tfsee.t <- t(tfsee)
tfsee.t <- as.data.frame(tfsee.t)
tfsee.t$rowMeans <- rowMeans(as.matrix(tfsee.t))
tfsee.t <- dplyr::filter(tfsee.t,rowMeans > summary(tfsee.t$rowMeans)[4])
tfsee.t <- tfsee.t[,-12]# Remove rowMeans column
tfsee <- t(tfsee.t)
write.csv(data.frame(colnames(tfsee)),"TFs_for_canSAR.csv")
tfsee <- t(scale(t(tfsee)))

# Go to https://cansarblack.icr.ac.uk/cpat and input list of TFs

# Fold in druggability information for TFs and read in canSAR data
drug.scores <- read.csv("./cpat-results.csv")
drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
colnames(drug.scores) <- c("TF","Score")
dim(drug.scores)

drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
dim(drug.scores)
drug.scores <- dplyr::arrange(drug.scores,TF)
drug.scores <- drug.scores[-9,]# Remove duplicate DDIT3
drug.scores <- drug.scores[drug.scores$TF %in% colnames(tfsee),]
drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)
length(which(drug.scores$TF== colnames(tfsee)))


# Rank frequency plot for subclones of endometrial metastasis
rownames(tfsee)
# Build TFSEE difference vector
tfsee.group1 <- tfsee[3,]
tfsee.group1 <- as.data.frame(tfsee.group1)
colnames(tfsee.group1) <- "group1"

tfsee.group2 <- tfsee[4,]
tfsee.group2 <- as.data.frame(tfsee.group2)
colnames(tfsee.group2) <- "group2"

print("Comparison: Subclone 1 v. Subclone 2")
print(rownames(tfsee)[3])
print(rownames(tfsee)[4])
print("Beginning comparison...")

length(which(rownames(tfsee.group1) == rownames(tfsee.group2)))
tfsee.group1 <- tfsee.group1
tfsee.group2 <- tfsee.group2
diff.tfsee <-tfsee.group1-tfsee.group2
hist(diff.tfsee)

colnames(diff.tfsee) <- "difference.tfsee"

total <- diff.tfsee

total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]

# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable

total <- dplyr::arrange(total,difference.tfsee)
total$rank <- 1:nrow(total)

p1<-ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_Sublcone_Comparison.pdf",width =5,height = 4)



# Rank frequency plot for carcinosarcoma
rownames(tfsee)
# Build TFSEE fold change vector
tfsee.group1 <- tfsee[7,]
tfsee.group1 <- as.data.frame(tfsee.group1)
colnames(tfsee.group1) <- "group1"

tfsee.group2 <- tfsee[8,]
tfsee.group2 <- as.data.frame(tfsee.group2)
colnames(tfsee.group2) <- "group2"


print("Comparison: Sarcoma v. Carcinoma")
print(rownames(tfsee)[7])
print(rownames(tfsee)[8])
print("Beginning comparison...")

length(which(rownames(tfsee.group1) == rownames(tfsee.group2)))
tfsee.group1 <- tfsee.group1
tfsee.group2 <- tfsee.group2
diff.tfsee <-tfsee.group1-tfsee.group2
hist(diff.tfsee)

colnames(diff.tfsee) <- "difference.tfsee"

total <- diff.tfsee

total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]


# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable

total <- dplyr::arrange(total,difference.tfsee)
total$rank <- 1:nrow(total)

p2<-ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_Carcinosarcoma.pdf",width =5,height = 4)




# Rank frequency plot for endometrial cancer of different histologies
rownames(tfsee)
# Build TFSEE fold change vector
tfsee.group1 <- colMeans(tfsee[3:4,])
tfsee.group1 <- as.data.frame(tfsee.group1)
colnames(tfsee.group1) <- "group1"

tfsee.group2 <- colMeans(tfsee[1:2,])
tfsee.group2 <- as.data.frame(tfsee.group2)
colnames(tfsee.group2) <- "group2"

print("Comparison: Serous v. Endometrioid EC")
print(paste0("Mean of:"))
print(rownames(tfsee[c(3,4),]))
print(paste0("Mean of:"))
print(rownames(tfsee[c(1,2),]))
print("Beginning comparison...")

length(which(rownames(tfsee.group1) == rownames(tfsee.group2)))
tfsee.group1 <- tfsee.group1
tfsee.group2 <- tfsee.group2
diff.tfsee <-tfsee.group1-tfsee.group2
hist(diff.tfsee)

colnames(diff.tfsee) <- "difference.tfsee"

total <- diff.tfsee

total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]


# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable

total <- dplyr::arrange(total,difference.tfsee)
total$rank <- 1:nrow(total)

p3<-ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_SerousVEndometrioid_Endometrial.pdf",width =5,height = 4)





# Rank frequency plot for ovarian cancer of different histologies
rownames(tfsee)
# Build TFSEE fold change vector
tfsee.group1 <- colMeans(tfsee[c(6,11),])
tfsee.group1 <- as.data.frame(tfsee.group1)
colnames(tfsee.group1) <- "group1"

tfsee.group2 <- tfsee[5,]
tfsee.group2 <- as.data.frame(tfsee.group2)
colnames(tfsee.group2) <- "group2"

print("Comparison: Serous v. Endometrioid OC")
print(paste0("Mean of:"))
print(rownames(tfsee[c(6,11),]))
print(rownames(tfsee)[5])
print("Beginning comparison...")

length(which(rownames(tfsee.group1) == rownames(tfsee.group2)))
tfsee.group1 <- tfsee.group1
tfsee.group2 <- tfsee.group2
diff.tfsee <-tfsee.group1-tfsee.group2
hist(diff.tfsee)

colnames(diff.tfsee) <- "difference.tfsee"

total <- diff.tfsee

total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]


# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable

total <- dplyr::arrange(total,difference.tfsee)
total$rank <- 1:nrow(total)

p4<-ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_SerousVEndometrioid_Ovarian.pdf",width =5,height = 4)



# Rank frequency plot for GIST versus HGSOC
rownames(tfsee)

# Build TFSEE fold change vector
tfsee.group1 <- colMeans(tfsee[c(9,10),])
tfsee.group1 <- as.data.frame(tfsee.group1)
colnames(tfsee.group1) <- "group1"

tfsee.group2 <- colMeans(tfsee[c(6,11),])
tfsee.group2 <- as.data.frame(tfsee.group2)
colnames(tfsee.group2) <- "group2"

print("Comparison: GIST v. HGSOC")
print(paste0("Mean of:"))
print(rownames(tfsee[c(9,10),]))
print(paste0("Mean of:"))
print(rownames(tfsee[c(6,11),]))
print("Beginning comparison...")

length(which(rownames(tfsee.group1) == rownames(tfsee.group2)))
tfsee.group1 <- tfsee.group1
tfsee.group2 <- tfsee.group2
diff.tfsee <-tfsee.group1-tfsee.group2
hist(diff.tfsee)

colnames(diff.tfsee) <- "difference.tfsee"

total <- diff.tfsee

total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]


# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable

total <- dplyr::arrange(total,difference.tfsee)
total$rank <- 1:nrow(total)

p5<-ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_GISTvHGSOC.pdf",width =5,height = 4)


CombinePlots(list(p4,p3,p5),ncol=1)+ggsave("Freq_Distribtions_Suppl.pdf",
                                           width=5,height = 12)

CombinePlots(list(p1,p2),ncol=1)+ggsave("Freq_Distribtions_Main.pdf",
                                           width=5,height =8)


# Plot heatmap with druggability
# Make heatmap annotation
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",5),rep("2",6))),
  Type = factor(c(rep("1",5),rep("2",5),rep("3",1))),
  Histology = factor(c(rep("1",5),rep("2",4),rep("3",1),rep("4",1))),
  Stage = factor(c(rep("1",1),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:11))))


library(ComplexHeatmap)

# Z-score the cell types
pdf(paste0("TFSEE_celltype_scaled.pdf"), width = 7, height =11)


drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)

tf.names <- colnames(tfsee)[colnames(tfsee) %in% drug.scores.sub$TF]
idx <- match(tf.names,colnames(tfsee))

ha.2 = rowAnnotation(foo = anno_mark(at = idx, labels = tf.names))

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

activity=list(Druggability=c("Yes"="Black","No"="gray60")) 
ha.3 = HeatmapAnnotation(
  Druggability = drug.scores$Druggable,which = "row",col =activity)


Heatmap(scale(t(tfsee)),top_annotation = ha, 
        show_row_names =T,clustering_method_rows = "complete",clustering_method_columns = "complete",clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        right_annotation = c(ha.2,ha.3))
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo_TFSEE_Plotting.txt")
```


Interested in more exciting research in cancer genomics? Visit https://www.thefrancolab.org/ to learn more!


