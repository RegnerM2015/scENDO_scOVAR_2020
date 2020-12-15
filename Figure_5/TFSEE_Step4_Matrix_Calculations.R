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
rna.pseudobulk <- rna.pseudobulk[enhancer.motifs$X,]
rna.pseudobulk[rna.pseudobulk == 0] <- NA
rna.pseudobulk <- na.omit(rna.pseudobulk)
dim(rna.pseudobulk)


TFs.to.use <- intersect(enhancer.motifs$X,rownames(rna.pseudobulk))
dim(enhancer.motifs)
enhancer.motifs <- enhancer.motifs[enhancer.motifs$X %in% TFs.to.use,]
dim(enhancer.motifs)


# Rename syntax of peaks
colnames(enhancer.motifs) <- sub("\\.",":",colnames(enhancer.motifs))
colnames(enhancer.motifs) <- sub("\\.","-",colnames(enhancer.motifs))

# Only include TFs in TF expression matrix
dim(enhancer.motifs)
enhancer.motifs <- enhancer.motifs[enhancer.motifs$X %in% rownames(rna.pseudobulk),]
dim(enhancer.motifs)

length(which(rownames(rna.pseudobulk) == enhancer.motifs$X))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.rna <- rlog(dds,blind = T)
saveRDS(rlog.rna,"rlog_rna.rds")



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

saveRDS(rlog.atac,"rlog_atac.rds")





###############################################
# Part 4: Scaling of each matrix from 0-1
###############################################

# Scale Enhancer matrix :
enhancer.mat <- assay(rlog.atac)
library(scales)
for (i in 1:ncol(enhancer.mat)){
  enhancer.mat[,i] <- rescale(enhancer.mat[,i],to = c(0,1))
}

#saveRDS(enhancer.motifs,"enhancer_motifs.rds")
# Scale TF motif matrix (invert the pvalue):
motif.mat <- enhancer.motifs
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

library(scales)
for (i in 1:ncol(expr.mat)){
  expr.mat[,i] <- rescale(expr.mat[,i],to = c(0,1))
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

pdf(paste0("TFSEE_TF_scaled_Zscore.pdf"), width = 5, height = 6)

Heatmap(t(scale(tfsee)),show_row_names =T)

dev.off()



# Z-score the cell types to reveal patterns/groups of TF/enhancers

pdf(paste0("TFSEE_celltype_scaled_Zscore.pdf"), width = 5, height = 6)

Heatmap(scale(t(tfsee)),show_row_names =T)

dev.off()


# Fold in druggability information for TFs
drug.scores <- read.csv("./CanSar_TF_data.csv")
drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
colnames(drug.scores) <- c("TF","Score")
dim(drug.scores)

drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
dim(drug.scores)
drug.scores <- dplyr::arrange(drug.scores,TF)
drug.scores <- drug.scores[-15,]# Remove duplicate DDIT3

length(which(drug.scores$TF== colnames(tfsee)))


# Pull out groups of TFs from heatmap based on drugability 
#####################################################
# Label the most druggable TFs

# Z-score the TFSEEs to reveal patterns/groups of cell types 

pdf(paste0("TFSEE_TF_scaled_Zscore_Drugs_labeled.pdf"), width = 8, height =10)

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

pdf(paste0("TFSEE_CellType_scaled_Zscore_Drugs_labeled.pdf"), width = 7, height =12)

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



# Scatter plot of log2FC TFSEE score versus log2FC TF expression (Serous Ovarian v. Endometrioid Ovarian)

# Subset data to the most robustly expressed TFs across pseudobulk clusters
expr <- as.data.frame(assay(rlog.rna))
expr$mean <- rowMeans(as.matrix(expr))
expr <- dplyr::arrange(expr,desc(mean))
expr <- dplyr::filter(expr,mean > summary(expr$mean)[2])
tfs <- rownames(expr)

tfsee <- tfsee[,colnames(tfsee) %in% tfs]

# Build TFSEE fold change vector
tfsee.serous <- t(tfsee[c(6,11),])
tfsee.serous <- as.data.frame(rowMeans(tfsee.serous))#Take mean TFSEE value between two serous clusters
colnames(tfsee.serous) <- "Serous"

tfsee.endo <- as.data.frame(tfsee[5,])
colnames(tfsee.endo) <- "Endometrioid"

length(which(rownames(tfsee.serous) == rownames(tfsee.endo)))

FC.tfsee <- log2(tfsee.serous$Serous/tfsee.endo$Endometrioid)
hist(FC.tfsee)

# Build Expression fold change vector 
expr <- assay(rlog.rna)
expr <- expr[rownames(expr) %in% tfs,]
expr.serous <- expr[,c(6,11)]
expr.serous <- as.data.frame(rowMeans(expr.serous))#Take mean TFSEE value between two serous clusters
colnames(expr.serous) <- "Serous"

expr.endo <- as.data.frame(expr[,5])
colnames(expr.endo) <- "Endometrioid"

length(which(rownames(expr.serous) == rownames(expr.endo)))

FC.expr <- log2(expr.serous$Serous/expr.endo$Endometrioid)
hist(FC.expr)

length(which(rownames(tfsee.serous) == rownames(expr.endo)))
length(which(rownames(tfsee.serous) == rownames(expr.serous)))
length(which(rownames(tfsee.endo) == rownames(expr.endo)))
length(which(rownames(tfsee.endo) == rownames(expr.serous)))


plot(FC.tfsee,FC.expr)


total <- as.data.frame(cbind(FC.expr,FC.tfsee))
rownames(total) <- rownames(tfsee.endo)
total$sum <- rowSums(total)


#total <- total[-27,] # Remove ESRRG as it had a zero TFSEE value causing division by zero



total.serous <- dplyr::filter(total,sum > 0)
total.endo <- dplyr::filter(total,sum < 0)

drug.serous <- intersect(rownames(total.serous),drug.scores.sub$TF)
drug.endo <- intersect(rownames(total.endo),drug.scores.sub$TF)

total.serous <- total.serous[rownames(total.serous) %in% drug.serous,]
total.endo <- total.endo[rownames(total.endo) %in% drug.endo,]



# Make Scatter plot for FC TFSEE versus FC Expression

# Add druggability to DF:
drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))
length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable
ggplot(total,aes(x=FC.tfsee,y = FC.expr))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("log2FC(TF Expression)")+theme_bw()+
  xlab("log2FC(Raw TFSEE Score)")+
  scale_color_manual(values = c("gray60","black"))+ggsave("TFSEE_FC.pdf",width =5,height = 4)



# # Intersect TFSEE enhancers with peak2gene links from full cohort analysis:
# 
# 
# # Can we use TFSEE to prioritze which Peak2Gene links we can followup on in BAM reads?
# 
# # Read in LAPTM4B peaks:
# test.peaks <- readRDS("./P2G_Hits_LAPTM4B_data.rds")
# test.peaks <- as.data.frame(test.peaks$peak)
# colnames(test.peaks) <- "peak"
# test.peaks <- tidyr::separate(test.peaks,peak,into = c("seqnames","start"),sep = ":" )
# test.peaks <- tidyr::separate(test.peaks,start,into = c("start","end"),sep = "-")
# ##################################################################################
# 
# tfsee.peaks <- as.data.frame(colnames(motif.mat))
# colnames(tfsee.peaks)<- "start"
# tfsee.peaks <- tidyr::separate(data = tfsee.peaks, col = start, into = c("start", "end"), sep = "\\-")
# tfsee.peaks <- tidyr::separate(data = tfsee.peaks, col = start, into = c("seqnames", "start"), sep = "\\:")
# 
# # Use bedtoolsr to intersect tfsee.peaks with P2G links for select gene
# intersect <- bedtoolsr::bt.intersect(test.peaks,tfsee.peaks,wa = T)
# colnames(intersect) <- c("seqnames","start","end")
# 
# intersect$coord <- paste0(intersect$seqnames,":",intersect$start,"-",intersect$end)
# 
# 
# 
# # Are any P2G links in tfsee matrix?
# sig.distal <- readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/TFSEE_Methods/Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")
# sig.distal <- dplyr::filter(sig.distal,EmpFDR < 0.05)
# sig.distal <- sig.distal[,3:5]
# colnames(sig.distal) <- c("seqnames","start","end")
# intersect <- bedtoolsr::bt.intersect(sig.distal,tfsee.peaks,wa = T)
# colnames(intersect) <- c("seqnames","start","end")
# 
# 
# 
# sig.distal <- readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/TFSEE_Methods/Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")
# sig.distal <- dplyr::filter(sig.distal,coord %in% paste0(intersect$seqnames,":",intersect$start,"-",intersect$end))
# 
# 
# p2g <- plotPeak2GeneHeatmap(atac,
#                             corCutOff = 0,
#                             FDRCutOff = 1,
#                             groupBy = "predictedGroup_ArchR",
#                             returnMatrices = T)# Remake P2G with matrices
# 
# mat <- p2g$RNA$matrix
# colnames(mat) <- make.unique(p2g$RNA$colData$groupBy)
# rownames(mat) <- make.unique(p2g$Peak2GeneLinks$gene)
# 
# kmeans <- p2g$RNA$kmeansId
# idx <- grep(i,rownames(mat))
# kmeans <- kmeans[idx]
# 
# peaks.genes <- p2g$Peak2GeneLinks
# peaks.genes <- peaks.genes[idx,]
# 
# peaks.genes <- as.data.frame(peaks.genes)
# 
# peaks.genes.sub <- dplyr::filter(peaks.genes,peak %in% colnames(motif.mat))

# 
# # Serous endometrial versus endometroid endometrial
# 
# # Build TFSEE fold change vector
# tfsee.serous <- t(tfsee[c(3:4),])
# tfsee.serous <- as.data.frame(rowMeans(tfsee.serous))#Take mean TFSEE value between two serous clusters
# colnames(tfsee.serous) <- "Serous"
# 
# tfsee.endo <- t(tfsee[1:2,])
# tfsee.endo <- as.data.frame(rowMeans(tfsee.endo))#Take mean TFSEE value between two serous clusters
# colnames(tfsee.endo) <- "Endometrioid"
# 
# length(which(rownames(tfsee.serous) == rownames(tfsee.endo)))
# 
# FC.tfsee <- log2(tfsee.serous$Serous/tfsee.endo$Endometrioid)
# hist(FC.tfsee)
# 
# # Build Expression fold change vector 
# expr <- assay(rlog.rna)
# expr.serous <- expr[,c(3:4)]
# expr.serous <- as.data.frame(rowMeans(expr.serous))#Take mean TFSEE value between two serous clusters
# colnames(expr.serous) <- "Serous"
# 
# expr.endo <- as.data.frame(expr[,1:2])
# expr.endo <- as.data.frame(rowMeans(expr.endo))#Take mean TFSEE value between two serous clusters
# colnames(expr.endo) <- "Endometrioid"
# 
# length(which(rownames(expr.serous) == rownames(expr.endo)))
# 
# FC.expr <- log2(expr.serous$Serous/expr.endo$Endometrioid)
# hist(FC.expr)
# 
# length(which(rownames(tfsee.serous) == rownames(expr.endo)))
# length(which(rownames(tfsee.serous) == rownames(expr.serous)))
# length(which(rownames(tfsee.endo) == rownames(expr.endo)))
# length(which(rownames(tfsee.endo) == rownames(expr.serous)))
# 
# 
# plot(FC.tfsee,FC.expr)
# 
# 
# total <- as.data.frame(cbind(FC.expr,FC.tfsee))
# rownames(total) <- rownames(tfsee.endo)
# total$sum <- rowSums(total)
# 
# 
# total <- total[-27,] # Remove ESRRG as it had a zero TFSEE value causing division by zero
# 
# 
# 
# total.serous <- dplyr::filter(total,sum < 0)
# total.endo <- dplyr::filter(total,sum > 0)
# 
# drug.serous <- intersect(rownames(total.serous),drug.scores.sub$TF)
# drug.endo <- intersect(rownames(total.endo),drug.scores.sub$TF)
# 
# total.serous <- total.serous[rownames(total.serous) %in% drug.serous,]
# total.endo <- total.endo[rownames(total.endo) %in% drug.endo,]
# 
# 
# 
# # Make Scatter plot for FC TFSEE versus FC Expression
# 
# # Add druggability to DF:
# drug.scores.new <- dplyr::filter(drug.scores,TF != "ESRRG")
# length(which(drug.scores.new$TF == rownames(total)))
# 
# total$drug.score <- drug.scores.new$Druggable
# ggplot(total,aes(x=FC.tfsee,y = FC.expr))+
#   geom_point(aes(color = drug.score),size = 2)+
#   ylab("log2FC(TF Expression)")+
#   scale_color_manual(values = c("gray60","black"))+
#   scale_x_continuous(name="log2FC(TFSEE)", 
#                      limits=c(-3, 6),
#                      breaks = c(-3,-2,-1,0,1,2,3,4,5,6))+theme_bw()+ggsave("TFSEE_FC.pdf",width = 6,height = 4)
# 
# 
# 
# 
# 
# 
# # Serous endometrial versus serous ovarian
# 
# # Build TFSEE fold change vector
# tfsee.serous <- t(tfsee[c(6,11),])
# tfsee.serous <- as.data.frame(rowMeans(tfsee.serous))#Take mean TFSEE value between two serous clusters
# colnames(tfsee.serous) <- "Serous"
# 
# tfsee.endo <- t(tfsee[3:4,])
# tfsee.endo <- as.data.frame(rowMeans(tfsee.endo))#Take mean TFSEE value between two serous clusters
# colnames(tfsee.endo) <- "Endometrioid"
# 
# length(which(rownames(tfsee.serous) == rownames(tfsee.endo)))
# 
# FC.tfsee <- log2(tfsee.endo$Endometrioid/tfsee.serous$Serous)
# hist(FC.tfsee)
# 
# # Build Expression fold change vector 
# expr <- assay(rlog.rna)
# expr.serous <- expr[,c(6,11)]
# expr.serous <- as.data.frame(rowMeans(expr.serous))#Take mean TFSEE value between two serous clusters
# colnames(expr.serous) <- "Serous"
# 
# expr.endo <- as.data.frame(expr[,3:4])
# expr.endo <- as.data.frame(rowMeans(expr.endo))#Take mean TFSEE value between two serous clusters
# colnames(expr.endo) <- "Endometrioid"
# 
# length(which(rownames(expr.serous) == rownames(expr.endo)))
# 
# FC.expr <- log2(expr.endo$Endometrioid/expr.serous$Serous)
# hist(FC.expr)
# 
# length(which(rownames(tfsee.serous) == rownames(expr.endo)))
# length(which(rownames(tfsee.serous) == rownames(expr.serous)))
# length(which(rownames(tfsee.endo) == rownames(expr.endo)))
# length(which(rownames(tfsee.endo) == rownames(expr.serous)))
# 
# 
# plot(FC.tfsee,FC.expr)
# 
# 
# total <- as.data.frame(cbind(FC.expr,FC.tfsee))
# rownames(total) <- rownames(tfsee.endo)
# total$sum <- rowSums(total)
# 
# 
# total <- total[-27,] # Remove ESRRG as it had a zero TFSEE value causing division by zero
# 
# 
# 
# total.serous <- dplyr::filter(total,sum < 0)
# total.endo <- dplyr::filter(total,sum > 0)
# 
# drug.serous <- intersect(rownames(total.serous),drug.scores.sub$TF)
# drug.endo <- intersect(rownames(total.endo),drug.scores.sub$TF)
# 
# total.serous <- total.serous[rownames(total.serous) %in% drug.serous,]
# total.endo <- total.endo[rownames(total.endo) %in% drug.endo,]
# 
# 
# 
# # Make Scatter plot for FC TFSEE versus FC Expression
# 
# # Add druggability to DF:
# drug.scores.new <- dplyr::filter(drug.scores,TF != "ESRRG")
# length(which(drug.scores.new$TF == rownames(total)))
# 
# total$drug.score <- drug.scores.new$Druggable
# ggplot(total,aes(x=FC.tfsee,y = FC.expr))+
#   geom_point(aes(color = drug.score),size = 2)+
#   ylab("log2FC(TF Expression)")+
#   scale_color_manual(values = c("gray60","black"))+
#   scale_x_continuous(name="log2FC(TFSEE)", 
#                      limits=c(-3, 6),
#                      breaks = c(-3,-2,-1,0,1,2,3,4,5,6))+theme_bw()+ggsave("TFSEE_FC.pdf",width = 6,height = 4)
# 
# 



