library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(ArchR)
library(stringr)
library(stringi)
library(dplyr)
library(tidyr)
h5disableFileLocking()

rna <- readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/Figure_Making/Figure_1/endo_ovar_All_scRNA_processed.rds")
# Plot cell type UMAPs for RNA/ATAC

# RNA
levels(factor(rna$RNA_snn_res.0.7))

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cell.type <- rna$cell.type
#Manually annotate 23-cluster as smooth muscle
rna.df$cell.type <- str_replace_all(rna.df$cell.type,"23-Stromal fibroblast","23-Smooth muscle cells")

rna.df$cluster <- as.factor(rna.df$cluster)
rna.df$cell.type <- as.factor(rna.df$cell.type)
levels(rna.df$cluster)
levels(rna.df$cell.type)


# Help sort the cluster numbers:
###############################
epi <- grep("pitheli",levels(rna.df$cell.type))
epi.new <- grep("-Ciliated",levels(rna.df$cell.type))
epi <- c(epi,epi.new)

fibro <- grep("ibro",levels(rna.df$cell.type))

smooth <- grep("mooth",levels(rna.df$cell.type))

endo <- grep("dothel",levels(rna.df$cell.type))

t.nk <- grep("T cell",levels(rna.df$cell.type))
t.nk.new <- grep("Lym",levels(rna.df$cell.type))
t.nk <- c(t.nk,t.nk.new)

mac <- grep("acrophage",levels(rna.df$cell.type))

mast <- grep("Mast",levels(rna.df$cell.type))

b <- grep("B cell",levels(rna.df$cell.type))


cell.types.idx <- c(epi,fibro,smooth,endo,t.nk,mac,mast,b)

store <- numeric(0)
for(i in 1:length(cell.types.idx)){
  name <- levels(rna.df$cell.type)[cell.types.idx[i]]
  print(gsub("-.*","",name))
  new.name <- gsub("-.*","",name)
  new.num <- as.numeric(new.name)
  store[i] <- new.num
}
print(store)
#####################################################

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
rna.df$cluster.new <- factor(x = rna.df$cluster, levels = my_levels)


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

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_RNA.pdf",width = 8,height = 7)

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_RNA-labels.pdf",width = 8,height = 7)



# Suppl. Figure CNV plot
# Plot CNV box: 
meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7,levels = rev(c(11,20,21,22,31,
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
                                                           28,35)))

ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot()+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  ggsave("CNV_BoxPlot.pdf",width = 4,height = 8)

# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "MUC16",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=MUC16))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CA125 expression")+
  scale_fill_manual(values=rev(cols))+
  coord_cartesian(xlim=c(0, 1.5))+
  ggsave("CA125_BoxPlot.pdf",width = 4,height = 8)
# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "WFDC2",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=WFDC2))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("HE4 expression")+
  scale_fill_manual(values=rev(cols))+
  coord_cartesian(xlim=c(0, 1.5))+
  ggsave("HE4_BoxPlot.pdf",width = 4,height = 8)
# Make violin plots for CD117 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "KIT",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=KIT))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CD117 expression")+
  scale_fill_manual(values=rev(cols))+
  ggsave("CD117_BoxPlot.pdf",width = 4,height = 8)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
