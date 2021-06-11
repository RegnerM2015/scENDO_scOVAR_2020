###########################################################
# Matt Regner
# Franco Lab
# Date: May 2020-May 2021
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) Find distal enhancer elements that have 
#            nonzero counts across all samples (clusters)
#         2) Subset Marker Enhancer BED files to new enhancers
###########################################################
library(ArchR)
library(stringr)
library(utils)
library(dplyr)
ArchR::addArchRThreads(threads = 32)
h5disableFileLocking()
enhancers <- read.delim("Marker_Enhancers_ArchR-nodups.bed",header = F)
enhancers <- paste0(enhancers$V1,":",enhancers$V2,"-",enhancers$V3)

# Pseudobulk ATAC enhancer matrix
proj.archr <- readRDS("./final_archr_proj_archrGS.rds")
peaks <- getPeakSet(proj.archr)

peak.mat <- getMatrixFromProject(proj.archr,useMatrix = "PeakMatrix")
peak.names <- paste0(peaks@seqnames,":",peaks@ranges)


cell.names <- colnames(assay(peak.mat))
peak.names <-peak.names


peak.mat <- assay(peak.mat)
colnames(peak.mat) <- cell.names
rownames(peak.mat) <- peak.names
head(peak.mat[,1:4])


# Only tumor cell clusters that are 100% patient specific 
df <- as.data.frame(proj.archr@cellColData) %>% dplyr::group_by(predictedGroup_ArchR) %>% 
  dplyr::count(Sample) 

idents <- table(df$predictedGroup_ArchR)
idents.sub <- idents[idents ==1]
idents.sub <- names(idents.sub)
idents.sub <- idents.sub[c(1,2,3,4,5,6,10,12,13,14,15)]
idents.sub
print(idents.sub)
##################################################################

peaks.pseudobulk <- data.frame(rownames(peak.mat))
for (i in idents.sub){
  cells <- rownames(dplyr::filter(as.data.frame(proj.archr@cellColData),predictedGroup_ArchR == i))
  
  peak.mat.sub <- peak.mat[,colnames(peak.mat) %in% cells]
  
  peaks.bulk <- Matrix::rowSums(peak.mat.sub)
  
  peaks.pseudobulk$i <- peaks.bulk
  colnames(peaks.pseudobulk)[dim(peaks.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}

rownames(peaks.pseudobulk) <- peaks.pseudobulk[,1]
peaks.pseudobulk <- peaks.pseudobulk[,-1]

peaks.pseudobulk <- peaks.pseudobulk[rownames(peaks.pseudobulk) %in% unique(enhancers),]

# Remove enhancer peaks that have ANY zero counts in ANY sample (no replicates)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]

saveRDS(peaks.pseudobulk,"Peak_matrix.rds")

writeLines(capture.output(sessionInfo()), "sessionInfo_TFSEE_Step1_NonZero_Enhancers.txt")