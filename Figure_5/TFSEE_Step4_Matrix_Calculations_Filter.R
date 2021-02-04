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


###############################################
# Part 1: Read in TF motif prediction matrix 
###############################################
# Read in TF motif prediction matrix
enhancer.motifs<- read.delim("./motif_enhancers.txt",sep = "\t")

# Set labels of cell types of interest

# Only use clusters with 100% patient specificity 
labels <- c("32-Unciliated epithelia 1",
                   "22-Unciliated epithelia 1",
                   "19-Epithelial cell",
                   "35-Epithelial cell",
                   "1-Epithelial cell",
                   "10-Epithelial cell",
                   "16-Fibroblast",
                   "18-Epithelial cell",
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


dim(rna.pseudobulk)
rna.pseudobulk.sub <- rna.pseudobulk[enhancer.motifs$X,]
rna.pseudobulk.sub[rna.pseudobulk.sub == 0] <- NA
rna.pseudobulk.sub <- na.omit(rna.pseudobulk.sub)
dim(rna.pseudobulk.sub)


TFs.to.use <- intersect(enhancer.motifs$X,rownames(rna.pseudobulk.sub))
dim(enhancer.motifs)
enhancer.motifs <- enhancer.motifs[enhancer.motifs$X %in% TFs.to.use,]
dim(enhancer.motifs)


# Rename syntax of peaks
colnames(enhancer.motifs) <- sub("\\.",":",colnames(enhancer.motifs))
colnames(enhancer.motifs) <- sub("\\.","-",colnames(enhancer.motifs))

# Only include TFs in TF expression matrix
dim(enhancer.motifs)
enhancer.motifs <- enhancer.motifs[enhancer.motifs$X %in% rownames(rna.pseudobulk.sub),]
dim(enhancer.motifs)

length(which(rownames(rna.pseudobulk.sub) == enhancer.motifs$X))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.rna <- rlog(dds,blind = T)
#saveRDS(rlog.rna,"rlog_rna.rds")



###############################################
# Part 3: Pseudobulk enhancer activity
###############################################

# Pseudobulk ATAC enhancer matrix
distal <- colnames(enhancer.motifs)[-1]

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

# Subset peak matrix to the enhancers of interest 
dim(peaks.mat)
peaks.mat <- peaks.mat[unique(distal),]
dim(peaks.mat)


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


dim(peaks.pseudobulk)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]
dim(peaks.pseudobulk)


rownames(enhancer.motifs) <- enhancer.motifs[,1]
enhancer.motifs <- enhancer.motifs[,-1]
head(enhancer.motifs[,1:4])
dim(enhancer.motifs)
enhancer.motifs <- enhancer.motifs[,colnames(enhancer.motifs) %in% rownames(peaks.pseudobulk)]
dim(enhancer.motifs)

length(rownames(peaks.pseudobulk) == colnames(enhancer.motifs))
# ######################################################

#DESeq ATAC

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = peaks.pseudobulk,
                              colData = data.frame(meta=colnames(peaks.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.atac <- rlog(dds,blind = T)

#saveRDS(rlog.atac,"rlog_atac.rds")





###############################################
# Part 4: Scaling of each matrix from 0-1
###############################################

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
library(scales)
for (i in 1:ncol(motif.mat)){
  motif.mat[,i] <- rescale(-log10(motif.mat[,i]+0.5),to = c(0,1))
}
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
  Patient = factor(as.character(seq(1:11)))
)


# Z-score the TFSEEs to reveal patterns/groups of cell types 

pdf(paste0("TFSEE_TF_scaled_Zscore_v2.pdf"), width = 5, height = 6)

Heatmap(t(scale(tfsee)),show_row_names =T)

dev.off()



# Z-score the cell types to reveal patterns/groups of TF/enhancers

pdf(paste0("TFSEE_celltype_scaled_Zscore_v2.pdf"), width = 5, height = 6)

Heatmap(scale(t(tfsee)),show_row_names =T)

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
