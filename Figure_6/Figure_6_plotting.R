
tfsee <- readRDS("./tfsee_matrix.rds")
tfsee.t <- t(tfsee)
tfsee.t <- as.data.frame(tfsee.t)
tfsee.t$rowMeans <- rowMeans(as.matrix(tfsee.t))
tfsee.t <- dplyr::filter(tfsee.t,rowMeans > summary(tfsee.t$rowMeans)[3])
tfsee.t <- tfsee.t[,-12]# Remove rowMeans column
tfsee <- t(tfsee.t)


tfsee <- t(scale(t(tfsee)))
expr.mat <- readRDS("./expr_matrix.rds")

expr.mat <- expr.mat[rownames(expr.mat) %in% colnames(tfsee),]
expr.mat <- expr.mat[ order(match(rownames(expr.mat), colnames(tfsee))), ]
length(which(rownames(expr.mat)==colnames(tfsee)))


# Fold in druggability information for TFs
drug.scores <- read.csv("./TFs_CanSAR.csv")
drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
colnames(drug.scores) <- c("TF","Score")
dim(drug.scores)

drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
dim(drug.scores)
drug.scores <- dplyr::arrange(drug.scores,TF)
drug.scores <- drug.scores[-26,]# Remove duplicate DDIT3
drug.scores <- drug.scores[drug.scores$TF %in% colnames(tfsee),]

drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)

length(which(drug.scores$TF== colnames(tfsee)))



# # Scatter plot of log2FC TFSEE score versus log2FC TF expression 

# Build TFSEE fold change vector
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

# Build Expression fold change vector

expr.group1 <- as.data.frame(expr.mat[,3])
colnames(expr.group1) <- "group1"
rownames(expr.group1) <- rownames(expr.mat)

expr.group2 <- as.data.frame(expr.mat[,4])
colnames(expr.group2) <- "group2"
rownames(expr.group2) <- rownames(expr.mat)

length(which(rownames(expr.group1) == rownames(expr.group2)))
expr.group1 <- expr.group1
expr.group2 <- expr.group2
diff.expr <- expr.group1-expr.group2
hist(diff.expr)

length(which(rownames(tfsee.group1) == rownames(expr.group2)))
length(which(rownames(tfsee.group1) == rownames(expr.group1)))
length(which(rownames(tfsee.group2) == rownames(expr.group2)))
length(which(rownames(tfsee.group2) == rownames(expr.group1)))


colnames(diff.tfsee) <- "difference.tfsee"
colnames(diff.expr) <- "difference.expr"
plot(diff.tfsee$difference.tfsee,diff.expr$difference.expr)


total <- as.data.frame(cbind(diff.expr,diff.tfsee))
total$sum <- rowSums(total)

total.group1 <- dplyr::filter(total,sum > 0)
total.group2 <- dplyr::filter(total,sum < 0)

drug.group1 <- intersect(rownames(total.group1),drug.scores.sub$TF)
drug.group2 <- intersect(rownames(total.group2),drug.scores.sub$TF)

total.group1 <- total.group1[rownames(total.group1) %in% drug.group1,]
total.group2 <- total.group2[rownames(total.group2) %in% drug.group2,]



# Make Scatter plot for FC TFSEE versus FC Expression

# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable
ggplot(total,aes(x=difference.tfsee,y = difference.expr))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("log2FC(TF Expression)")+theme_bw()+
  xlab("log2FC(Raw TFSEE Score)")+
  scale_color_manual(values = c("gray60","black"))+
  geom_vline(aes(xintercept = 0.5),linetype = "dashed") + 
  geom_vline(aes(xintercept = -0.5),linetype = "dashed") + 
  geom_hline(aes(yintercept = 0.5),linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.5),linetype = "dashed") + ggsave("TFSEE_FC_Scar_Epi_Ov.pdf",width =5,height = 4)

