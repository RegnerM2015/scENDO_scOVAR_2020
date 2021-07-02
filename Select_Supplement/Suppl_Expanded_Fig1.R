library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(stringr)
library(ArchR)


atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)
atac.df$cell.type <- as.factor(atac.df$cell.type)
atac.df$cell.type <- str_replace_all(atac.df$cell.type,"23-Stromal fibroblast","23-Smooth muscle cells")
atac.df$cluster <-  gsub("-.*","",atac.df$cell.type)


my_levels <- c(11,20,21,22,31,
               19,34,
               3,
               9,10,
               16,17,
               0,27,
               6,8,12,14,15,18,24,25,26,29,
               7,23,
               1,33,
               2,4,30,
               5,13,
               32,
               28,35)


# Relevel object@ident
atac.df$cluster.new <- factor(x = atac.df$cluster, levels = my_levels)


epithelial.cols <- colorRampPalette(c("#A0E989", "#245719"))
epithelial.cols <- epithelial.cols(14)

fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(10)

smooth.cols <- c("#b47fe5","#d8b7e8")

endo.cols <- c("#93CEFF","#4A99FF")

t.cols <- c("gray60","gray20","gray40")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)

p1 <- ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
p1 <- LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)

# Make pred UMAP
atac$barcode <- rownames(atac@cellColData)
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac$idx <- 1:nrow(atac@cellColData)

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

length(which(atac$predictedGroup_ArchR==atac.df$cluster))
atac.df$predictedScore_ArchR <- as.numeric(atac$predictedScore_ArchR)

p2 <- ggplot(atac.df,aes(x = x,y=y,color =predictedScore_ArchR))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_gradient2(midpoint=0.6,mid="white",low="blue",high="red3")
p2+ggsave("FullCohort_ATAC_predictionscore.pdf",width = 6,height = 4)


ggplot(atac.df,aes(x=predictedScore_ArchR))+
         geom_histogram(binwidth = 0.025,fill="gray", color="black", alpha=0.9)+
  theme_classic()+
  geom_vline(xintercept = 0.50)+ggsave("FullCohort_ATAC_predictionscore_histogram.pdf",width = 4,height = 4)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
