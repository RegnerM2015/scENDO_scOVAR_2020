library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(stringr)
library(stringi)
library(dplyr)
library(tidyr)

rna <- readRDS("endo_EEC_scRNA_processed.rds")
# Plot cell type UMAPs for RNA/ATAC

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type
rna.df$cell.type <- str_replace(rna.df$cell.type,pattern = "12-Stromal fibroblasts","12-Smooth muscle cells")
#Rename 12 cluster as smooth muscle based on Myh11 expression
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cluster.new <- factor(rna.df$cluster,levels = c("6",
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
))


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

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color =cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))+ggsave("Cell_type_RNA.pdf",width = 8,height = 6)
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
               "27","28")
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

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("6",
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
                                                      "27","28"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[1:5])+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(position = "right")+
  scale_x_discrete(position = "top")+NoLegend()+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells") -> pfirst



CombinePlots(list(p3,p2,pfirst),ncol=3)+ggsave("Malignant_markers_vln.pdf",width=12,height =6)

p.CNV <- FeaturePlot(rna,features = "Total_CNVs")
p.CA125 <- FeaturePlot(rna,features = "MUC16")
p.HE4 <- FeaturePlot(rna,features = "WFDC2")

CombinePlots(list(plabels,p.CA125,p.HE4),ncol=1)+ggsave("Malignant_Features.pdf",width = 6,height =16)



malignant <- c("6",
               "8",
               "10",
               "14",
               "15",
               "19",
               "21",
               "22")
remaining <- c("1",
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
               "27","28")
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
