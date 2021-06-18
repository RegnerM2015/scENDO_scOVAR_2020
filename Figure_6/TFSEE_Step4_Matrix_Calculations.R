###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: All 
# Description: This script performs the following tasks  
#     1) Read in TF motif prediction matrix 
#     2) Pseudobulk TF expression (DESeq rlog transformation)
#     3) Pseudobulk enhancer activity (DESeq rlog transformation)
#     4) Scaling of each matrix from 0-1
#     5) TFSEE matrix operations 
#################################################################
library(ArchR)

library(DESeq2)

###############################################
# Part 1: Read in TF motif prediction matrix 
###############################################
# Read in TF motif prediction matrix
enhancer.motifs<- read.delim("./motif_enhancers_scaled.txt",sep = "\t")
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
#saveRDS(rlog.rna,"rlog_rna.rds")


###############################################
# Part 3: Pseudobulk enhancer activity
###############################################

# Pseudobulk ATAC enhancer matrix

atac <- readRDS("./final_archr_proj_archrGS.rds")
peaks.mat <- getMatrixFromProject(atac,useMatrix = "PeakMatrix")
peaks <- getPeakSet(atac)

peak.names <- paste0(peaks@seqnames,":",peaks@ranges)

length(peak.names)
dim(peaks.mat)

cell.names <- colnames(assay(peaks.mat))
peak.names <- peak.names

peaks.mat <- assay(peaks.mat)
colnames(peaks.mat) <- cell.names
rownames(peaks.mat) <- peak.names
head(peaks.mat[,1:4])

# Construct pseudobulk peak/enhancer matrix
peaks.pseudobulk <- data.frame(rownames(peaks.mat))
for (i in labels){
  cells <- rownames(dplyr::filter(as.data.frame(atac@cellColData),predictedGroup_ArchR == i))
  
  peaks.sub <- peaks.mat[,colnames(peaks.mat) %in% cells]
  
  peaks.bulk <- rowSums(as.matrix(peaks.sub))
  
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

#saveRDS(rlog.atac,"rlog_atac.rds")



###############################################
# Part 4: Subsetting and scaling from 0 to 1
###############################################

# Subset and reorder motif matrix, peak matrix and expression matrix



# Scale Enhancer matrix :
enhancer.mat <- assay(rlog.atac)
enhancer.mat <- as.data.frame(enhancer.mat)
enhancer.mat$mean <- rowMeans(enhancer.mat)
enhancer.mat <- dplyr::filter(enhancer.mat,mean > summary(enhancer.mat$mean)[2])
enhancer.mat <- enhancer.mat[,-12]
enhancer.mat <- as.matrix(enhancer.mat)
library(scales)
for (i in 1:ncol(enhancer.mat)){
  enhancer.mat[,i] <- rescale(enhancer.mat[,i],to = c(0,1))
}

#saveRDS(enhancer.motifs,"enhancer_motifs.rds")
# Scale TF motif matrix (invert the pvalue):
motif.mat <- enhancer.motifs
motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(enhancer.mat)]
length(which(colnames(motif.mat) == rownames(enhancer.mat)))
# library(scales)
# for (i in 1:ncol(motif.mat)){
#   motif.mat[,i] <- rescale(-log10(motif.mat[,i]+0.5),to = c(0,1))
# }
# library(scales)
# for (i in 1:ncol(motif.mat)){
#   motif.mat[,i] <- rescale(motif.mat[,i],to = c(0,1))
# }


# Scale TF expression matrix:
expr.mat <- assay(rlog.rna)
expr.mat <- as.data.frame(expr.mat)
expr.mat <- expr.mat[rownames(expr.mat) %in% TFs.to.use,]
expr.mat$mean <- rowMeans(expr.mat)
expr.mat <- dplyr::filter(expr.mat,mean > summary(expr.mat$mean)[2])
expr.mat <- expr.mat[,-12]
expr.mat <- as.matrix(expr.mat)
library(scales)
for (i in 1:ncol(expr.mat)){
  expr.mat[,i] <- rescale(expr.mat[,i],to = c(0,1))
}

TFs <- intersect(rownames(motif.mat),rownames(expr.mat))
expr.mat <- expr.mat[TFs,]
motif.mat <- motif.mat[TFs,]
length(which(rownames(motif.mat) == rownames(expr.mat)))

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
saveRDS(expr.mat,"expr_matrix.rds")
saveRDS(motif.mat,"motif_matrix.rds")

tfsee <- readRDS("tfsee_matrix.rds")
##################################################################
# END OF TFSEE COMPUTATION
##################################################################



##################################################################
# Plotting of TFSEE heatmap and log2FC scatter plot 
##################################################################

# Make heatmap annotation
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",5),rep("2",6))),
  Type = factor(c(rep("1",5),rep("2",5),rep("3",1))),
  Histology = factor(c(rep("1",5),rep("2",4),rep("3",1),rep("4",1))),
  Stage = factor(c(rep("1",1),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:11))),
  size = DESeq2::sizeFactors(dds))


# Z-score the TFSEEs to reveal patterns/groups of cell types 

pdf(paste0("TFSEE_TF_scaled.pdf"), width = 5, height = 6)

Heatmap(t(scale(tfsee)),show_row_names =T,clustering_method_columns = "ward.D",top_annotation = ha)

dev.off()


pdf(paste0("TFSEE_TF_unscaled.pdf"), width = 5, height = 6)

Heatmap(t(tfsee),show_row_names =T,clustering_method_columns = "ward.D",top_annotation = ha,
        clustering_method_rows = "ward.D")

dev.off()

# Z-score the cell types to reveal patterns/groups of TF/enhancers

pdf(paste0("TFSEE_celltype_scaled.pdf"), width = 5, height = 9)

Heatmap(scale(t(tfsee)),top_annotation = ha,
        show_row_names =T,clustering_method_rows = "ward.D",clustering_method_columns = "ward.D")

dev.off()
pdf(paste0("TFSEE_celltype_scaled.pdf"), width = 5, height = 6)

Heatmap(scale(t(tfsee)),show_row_names =T,cluster_rows = T,cluster_columns = T,clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2")

dev.off()

# Fold in druggability information for TFs
drug.scores <- read.csv("./TFs_CanSAR.csv")
drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
colnames(drug.scores) <- c("TF","Score")
dim(drug.scores)

drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
dim(drug.scores)
drug.scores <- dplyr::arrange(drug.scores,TF)
drug.scores <- drug.scores[-18,]# Remove duplicate DDIT3
drug.scores <- drug.scores[drug.scores$TF %in% colnames(tfsee),]


length(which(drug.scores$TF== colnames(tfsee)))


# Pull out groups of TFs from heatmap based on drugability 
#####################################################
# Label the most druggable TFs

# Z-score the TFSEEs to reveal patterns/groups of cell types 

pdf(paste0("TFSEE_TF_scaled_Zscore_Drugs_labeled_v2.pdf"), width = 8, height =10)

drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)

tf.names <- colnames(tfsee)[colnames(tfsee) %in% drug.scores.sub$TF]
idx <- match(tf.names,colnames(tfsee))

ha.2 = rowAnnotation(foo = anno_mark(at = idx, labels = tf.names))

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, TRUE, FALSE)

ha.3 = rowAnnotation(
  Druggability = drug.scores$Druggable)

Heatmap(t(scale(tfsee)),
        top_annotation = ha,show_row_names =T,row_split = 4,right_annotation = c(ha.2,ha.3))

dev.off()

# Z-score the cell types to reveal patterns/groups of TFSEEs

pdf(paste0("TFSEE_CellType_scaled_Zscore_Drugs_labeled_v2.pdf"), width = 7, height =12)

drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)

tf.names <- colnames(tfsee)[colnames(tfsee) %in% drug.scores.sub$TF]
idx <- match(tf.names,colnames(tfsee))

ha.2 = rowAnnotation(foo = anno_mark(at = idx, labels = tf.names))

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")


activity=list(Druggability=c("Yes"="Black","No"="gray60")) 
ha.3 = HeatmapAnnotation(
  Druggability = drug.scores$Druggable,which = "row",col =activity)

Heatmap(scale(t(tfsee)),
        top_annotation = ha,show_row_names =T,right_annotation = c(ha.2,ha.3))

dev.off()


##########################
# END OF SCRIPT
##########################
