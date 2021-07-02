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
atac.df$cell.type <- str_replace_all(atac.df$cell.type,"12-Stromal fibroblast","12-Smooth muscle cells")
atac.df$cluster <-  gsub("-.*","",atac.df$cell.type)


my_levels <- c("6",
               "8",
               "10",
               "14",
               "15",
               "19",
               "21",
               "22",
               "1",
               "2",
               "9",
               "13",
               "16",
               "17",
               "20","23",
               "4" ,"5","12",
               "3","11","26",
               "0","18",
               "24","7",
               "25",
               "27","28"
)


# Relevel object@ident
atac.df$cluster.new <- factor(x = atac.df$cluster, levels = my_levels)




epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#f593f1","#f272ed","#ef52e9","#e415dd","#c412bd","#a30f9d","#820c7e","#62095e")
smooth.cols <- c("#c2abd6","#b47fe5","#8948a1")

endo.cols <- c("#93CEFF","#4A99FF","#286ab5")

t.cols <- c("gray40","gray60")

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
p2+ggsave("EEC_ATAC_predictionscore.pdf",width = 6,height = 4)


ggplot(atac.df,aes(x=predictedScore_ArchR))+
         geom_histogram(binwidth = 0.025,fill="gray", color="black", alpha=0.9)+
  theme_classic()+
  geom_vline(xintercept = 0.50)+ggsave("EEC_ATAC_predictionscore_histogram.pdf",width = 4,height = 4)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
