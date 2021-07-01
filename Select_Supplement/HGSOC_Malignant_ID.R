library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(stringr)
library(stringi)
library(dplyr)
library(tidyr)

rna <- readRDS("ovar_HGSOC_scRNA_processed.rds")
# Plot cell type UMAPs for RNA/ATAC

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cluster.new <- factor(rna.df$cluster,levels = c("0","2","3","7","11","15","16","19",
                                                       "1","4","9","10","22",
                                                       "14",
                                                       "5","12","18","21",
                                                       "6","8","13","20","23",
                                                       "17"
))


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281","#6B0772")

endo.cols <- "#377EB8"

t.cols <- c("#CCCCCC","#333333","#666666","#999999")

macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(5)


b.cols <- "#E41A1C"


cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color =cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))
rna.sample.plot +ggsave("Cell_type_RNA.pdf",width = 8,height = 6)
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



my_levels <- c("0","2","3","7","11","15","16","19",
               "1","4","9","10","22",
               "14",
               "5","12","18","21",
               "6","8","13","20","23",
               "17")
# Plot CNV box: 
meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7,levels = rev(c("0","2","3","7","11","15","16","19",
                                                           "1","4","9","10","22",
                                                           "14",
                                                           "5","12","18","21",
                                                           "6","8","13","20","23",
                                                           "17")))

p0<-ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot(lwd=0.45,outlier.size = 0.55,fatten=0.95)+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  scale_y_continuous(position = "right")
# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "MUC16",pt.size = 0)+coord_flip()+NoLegend()
p2 <- ggplot(p1$data,aes(y=ident,x=MUC16))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CA125 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "WFDC2",pt.size = 0)+coord_flip()+NoLegend()
p3 <- ggplot(p1$data,aes(y=ident,x=WFDC2))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("HE4 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")



library(forcats)

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])


library(dplyr)
# Patient proportion per subcluster in RNA:
###########################################################
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("0","2","3","7","11","15","16","19",
                                                      "1","4","9","10","22",
                                                      "14",
                                                      "5","12","18","21",
                                                      "6","8","13","20","23",
                                                      "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[c(8,11)])+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(position = "right")+
  scale_x_discrete(position = "top")+NoLegend()+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")-> pfirst




CombinePlots(list(p3,p2,p0,pfirst),ncol=4)+ggsave("Malignant_markers_vln.pdf",width=14,height =6)

p.CNV <- FeaturePlot(rna,features = "Total_CNVs")
p.CA125 <- FeaturePlot(rna,features = "MUC16")
p.HE4 <- FeaturePlot(rna,features = "WFDC2")

CombinePlots(list(plabels,p.CNV,p.CA125,p.HE4),ncol=1)+ggsave("Malignant_Features.pdf",width = 6,height =18)



malignant <- c("0","2","3","7","11","15","16","19")
remaining <- c("1","4","9","10","22",
               "14",
               "5","12","18","21",
               "6","8","13","20","23",
               "17")
# find markers avg logFC >= 0.5 and adj pval <=0.01
for ( i in malignant){
  
  tryCatch(
    markers <- FindMarkers(rna,ident.1=i,
                           ident.2=remaining,
                           features=c("MUC16",
                                      "WFDC2"),
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
