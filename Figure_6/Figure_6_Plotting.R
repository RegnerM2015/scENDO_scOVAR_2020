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

ggplot(total,aes(x=rank,y = difference.tfsee))+
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

ggplot(total,aes(x=rank,y = difference.tfsee))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("difference(Scaled TFSEE score)")+theme_bw()+
  xlab("Rank")+
  scale_color_manual(values = c("gray60","black"))+
  geom_hline(yintercept = 0.25,linetype="dashed")+
  geom_hline(yintercept = -0.25,linetype="dashed")+
  ggsave("TFSEE_Carcinosarcoma.pdf",width =5,height = 4)







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
