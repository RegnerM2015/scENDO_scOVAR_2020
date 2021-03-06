
###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) Find distal enhancer elements that have 
#            nonzero counts across all samples (clusters)
#         2) Subset Marker Enhancer BED files to new enhancers
###########################################################


# ONlY USE 100% patient-specfic 
labels <- list.files(pattern = "Marker_Enhancers_ArchR")
labels <- str_remove(labels,"Marker_Enhancers_ArchR_")
labels <- str_remove(labels,".bed")
labels <- labels[-c(12,13)]


# Read in input enhancers (PRIOR TO THIS R SCRIPT: bash cat command concatenated all files into one)
# cat Marker_Enhancers_ArchR* > Marker_Enhancers_ArchR.bed
# sort Marker_Enhancers_ArchR.bed | uniq -u > Marker_Enhancers_ArchR-nodups.bed

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

peaks.pseudobulk <- peaks.pseudobulk[rownames(peaks.pseudobulk) %in% unique(enhancers),]

# Remove enhancer peaks that have ANY zero counts in ANY sample (no replicates)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]



# Subset original Marker enhancer lists to new updated nonzero enhancers
for ( i in 1:length(labels)){
  file <- paste0("Marker_Enhancers_ArchR_",labels[i],".bed")
  
  bed <- read.delim(file,header = F)
  rownames(bed) <- paste0(bed$V1,":",bed$V2,"-",bed$V3)
  
  bed.new <- bed[rownames(bed) %in% rownames(peaks.pseudobulk),]
  
  write.table(bed.new,paste0("Marker_Enhancers_ArchR_",labels[i],"-updated.bed"),row.names = F,col.names = F,quote = F,sep = "\t")
}


