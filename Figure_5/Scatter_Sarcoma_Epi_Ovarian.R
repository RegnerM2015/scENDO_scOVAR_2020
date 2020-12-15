
tfsee <- readRDS("./tfsee_matrix.rds")
expr.mat <- readRDS("./expr_matrix.rds")

# Fold in druggability information for TFs
drug.scores <- read.csv("./CanSar_TF_data.csv")
drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
colnames(drug.scores) <- c("TF","Score")
dim(drug.scores)

drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
dim(drug.scores)
drug.scores <- dplyr::arrange(drug.scores,TF)
drug.scores <- drug.scores[-15,]# Remove duplicate DDIT3
drug.scores <- drug.scores[drug.scores$TF %in% colnames(tfsee),]

drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)

length(which(drug.scores$TF== colnames(tfsee)))



# # Scatter plot of log2FC TFSEE score versus log2FC TF expression (Sarcoma v. Epithelial of carsarc)

# Subset data to the most robustly expressed TFs across pseudobulk clusters
tfs <- rownames(expr.mat)

tfsee <- tfsee[,colnames(tfsee) %in% tfs]

# Build TFSEE fold change vector
tfsee.serous <- tfsee[7,]
tfsee.serous <- as.data.frame(tfsee.serous)
colnames(tfsee.serous) <- "Serous"

tfsee.endo <- tfsee[8,]
tfsee.endo <- as.data.frame(tfsee.endo)
colnames(tfsee.endo) <- "Endometrioid"

length(which(rownames(tfsee.serous) == rownames(tfsee.endo)))
# tfsee.serous <- tfsee.serous+0.1
# tfsee.endo <- tfsee.endo+0.1
FC.tfsee <- log2(tfsee.serous$Serous/tfsee.endo$Endometrioid)

hist(FC.tfsee)

# Build Expression fold change vector

expr.serous <- as.data.frame(expr.mat[,7])
colnames(expr.serous) <- "Serous"

expr.endo <- as.data.frame(expr.mat[,8])
colnames(expr.endo) <- "Endometrioid"

length(which(rownames(expr.serous) == rownames(expr.endo)))
# expr.serous <- expr.serous+0.1
# expr.endo <- expr.endo+0.1
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

total.serous <- dplyr::filter(total,sum > 0)
total.endo <- dplyr::filter(total,sum < 0)

drug.serous <- intersect(rownames(total.serous),drug.scores.sub$TF)
drug.endo <- intersect(rownames(total.endo),drug.scores.sub$TF)

total.serous <- total.serous[rownames(total.serous) %in% drug.serous,]
total.endo <- total.endo[rownames(total.endo) %in% drug.endo,]



# Make Scatter plot for FC TFSEE versus FC Expression

# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable
ggplot(total,aes(x=FC.tfsee,y = FC.expr))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("log2FC(TF Expression)")+theme_bw()+
  xlab("log2FC(Raw TFSEE Score)")+
  scale_color_manual(values = c("gray60","black"))+ggsave("TFSEE_FC_Scar_Epi_Ov.pdf",width =5,height = 4)

