library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(stringr)
library(stringi)
library(dplyr)
library(tidyr)

rna <- readRDS("endo_ovar_All_scRNA_processed.rds")
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
plabels <- LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_RNA-labels.pdf",width = 8,height = 7)


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

p0<-ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot(lwd=0.45,outlier.size = 0.55,fatten=0.95)+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  scale_y_continuous(position = "right")+
  ggsave("CNV_BoxPlot.pdf",width = 4,height = 8)

# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "MUC16",pt.size = 0)+coord_flip()+NoLegend()
p2 <- ggplot(p1$data,aes(y=ident,x=MUC16))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CA125 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")+
  ggsave("CA125_BoxPlot.pdf",width = 4,height = 8)
# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "WFDC2",pt.size = 0)+coord_flip()+NoLegend()
p3 <- ggplot(p1$data,aes(y=ident,x=WFDC2))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("HE4 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")+
  ggsave("HE4_BoxPlot.pdf",width = 4,height = 8)
# Make violin plots for CD117 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "KIT",pt.size = 0)+coord_flip()+NoLegend()
p4 <- ggplot(p1$data,aes(y=ident,x=KIT))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CD117 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")+
  ggsave("CD117_BoxPlot.pdf",width = 4,height = 8)


library(forcats)
sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])

# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c(11,20,21,22,31,
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
                                                      28,35))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+NoLegend()+
  scale_y_continuous(position = "right")+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 8) -> pfirst



CombinePlots(list(p4,p3,p2,p0,pfirst),ncol=5)+ggsave("Malignant_markers_vln.pdf",width=16,height =6)

p.CNV <- FeaturePlot(rna,features = "Total_CNVs")
p.CA125 <- FeaturePlot(rna,features = "MUC16")
p.HE4 <- FeaturePlot(rna,features = "WFDC2")
p.CD117 <- FeaturePlot(rna,features = "KIT")

CombinePlots(list(plabels,p.CNV,p.CA125,p.HE4,p.CD117),ncol=1)+ggsave("Malignant_Features.pdf",width = 6,height = 24)



malignant <- c(11,20,21,22,31,
               19,34,
               3,
               9,10,
               16,17,
               0,27)
remaining <- c(6,8,12,14,15,18,24,25,26,29,
               7,23,
               1,33,
               2,4,30,
               5,13,
               32,
               28,35)
# find markers avg logFC >= 0.5 and adj pval <=0.01
for ( i in malignant){
  
  tryCatch(
    markers <- FindMarkers(rna,ident.1=i,
                           ident.2=remaining,
                           features=c("MUC16",
                                      "WFDC2",
                                      "KIT"),
                           min.pct = 0,only.pos = T,
                           logfc.threshold = 0)
  )
  
  if(nrow(markers) >0){
    print(paste0("Cancer markers present for ",i,":"))
    print(markers)
    print(paste0("End markers for ",i))
  }else{
    
  }
  
}





writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
