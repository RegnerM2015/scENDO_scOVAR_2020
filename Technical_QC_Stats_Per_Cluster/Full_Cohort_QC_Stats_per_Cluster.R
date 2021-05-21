
###############################################################
# Matt Regner
# Franco Lab
#####################################################################

library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(stringr)
library(ArchR)

###########################################################
# Part 1: Original (no batch correction)
###########################################################

rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")
# Make patient sample metadata and color assignments 

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])



sampleAnnot <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L",
                                     "366C5L","37EACL","38FE7L","3BAE2L","3CCF1L","3E4D1L","3E5CFL"),
                          Color = sampleColors,
                          Cancer = c("endometrial","endometrial","endometrial","endometrial","endometrial","endometrial",
                                     "ovarian","ovarian","ovarian","ovarian","ovarian"),
                          Histology = c("endometrioid","endometrioid","endometrioid","endometrioid","endometrioid",
                                        "serous","endometrioid","serous","carcinosarcoma","GIST","serous"),
                          BMI = c(39.89,30.5,38.55,55.29,49.44,29.94,34.8,22.13,23.72,33.96,22.37),
                          Age = c(70,70,70,49,62,74,76,61,69,59,59),
                          Race = c("AA","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","AS"),
                          Stage = c("IA","IA","IA","IA","IA","IIIA","IA","IIB","IVB","IV","IIIC"),
                          Site = c("Endometrium","Endometrium","Endometrium","Endometrium","Endometrium","Ovary","Ovary","Ovary","Ovary","Ovary","Ovary"),
                          Type = c("Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Ovarian","Ovarian","Ovarian","Gastric","Ovarian"))



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
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
p1 <- LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)

rna$log10nCountRNA <- log10(rna$nCount_RNA)
p2 <- FeaturePlot(rna,features = "log10nCountRNA")
p3 <- FeaturePlot(rna,features = "nFeature_RNA")
CombinePlots(list(p1,p2,p3),ncol=3)+ggsave("Technical_Clusters_UMAP_RNA.pdf",width = 18,height = 6)
###########################################################################


# Plot nCountRNA and nFeatureRNA per tumor cluster
meta <- rna@meta.data
meta <- meta[meta$RNA_snn_res.0.7 %in% c(11,20,21,22,31,19,34,3,9,10,16,17,0,27),]

cells <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(meta$RNA_snn_res.0.7))))
names(cluster.ids) <- levels(factor(meta$RNA_snn_res.0.7))

for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)

levels(factor(meta$seurat_clusters))

cols <- test$Color

p1 <- ggplot(meta,aes(x=seurat_clusters,y=log10nCountRNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p2 <- ggplot(meta,aes(x=seurat_clusters,y=nFeature_RNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 500,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per epithelial cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scRNA_PercentMito_per_Cluster-epi.pdf",width = 8,height = 3)



# Plot nCountRNA and nFeatureRNA per fibroblast cluster
meta <- rna@meta.data
idx <- levels(factor(meta$cell.type))
idx <- grep("ibroblast",levels(factor(meta$cell.type)))
fibroblast <- levels(factor(meta$cell.type))[idx]
fibroblast <- fibroblast[-c(1,5,7,11)]# Remove non-normal fibroblasts
meta <- meta[meta$cell.type %in% fibroblast,]

cells <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(meta$RNA_snn_res.0.7))))
names(cluster.ids) <- levels(factor(meta$RNA_snn_res.0.7))

for ( i in factor(cells$RNA_snn_res.0.7)){
  
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)

levels(factor(meta$seurat_clusters))

cols <- test$Color



p3 <- ggplot(meta,aes(x=seurat_clusters,y=log10nCountRNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p4 <- ggplot(meta,aes(x=seurat_clusters,y=nFeature_RNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 500,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per fibroblast cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scRNA_PercentMito_per_Cluster-fibro.pdf",width = 8,height = 3)

CombinePlots(list(p1,p2,p3,p4),ncol=1)+ggsave("BoxPlots_scRNA_variables.pdf",width = 6,height = 12)


# Patient Proportion per Cluster
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
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 8)


###########################################################
# Part 2: RNA clusters with regressions
###########################################################

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
rna <- ScaleData(rna,features = VariableFeatures(rna),vars.to.regress = c("nCount_RNA","nFeature_RNA"))
rna <- RunPCA(rna)
rna <- FindNeighbors(rna,dims = 1:50)
rna <- FindClusters(rna,resolution = 0.7)
rna <- RunUMAP(rna,dims = 1:50)

# PanglaoDB 
tsv=gzfile("PanglaoDB_markers_27_Mar_2020.tsv.gz")  
panglaodb <- read.csv(tsv,header=T,sep = "\t") 
panglaodb <- dplyr::filter(panglaodb,species == "Hs" | species == "Mm Hs")# Human subset 
panglaodb <- dplyr::filter(panglaodb,organ == "Connective tissue" |
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle"
)
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type)

ESTIMATE.signatures <- "ESTIMATE_signatures.csv"

SAMPLE.ID = "endo_ovar_All"

# Verify SingleR annotations and check for mast cells:
rna <- AddModuleScore(rna,features = list(panglaodb$`B cells`,
                                          panglaodb$`Plasma cells`,
                                          panglaodb$`Mast cells`,
                                          panglaodb$Macrophages,
                                          panglaodb$`Dendritic cells`,
                                          panglaodb$`T cells`,
                                          panglaodb$`NK cells`,
                                          panglaodb$`Endothelial cells`,
                                          panglaodb$Fibroblasts,
                                          panglaodb$`Epithelial cells`,
                                          panglaodb$`Smooth muscle cells`,
                                          c("TPSB2","TPSAB1","KIT")),#Three gene Mast signature
                      name = c("B.","Plasma.","Mast.","Macrophage.","DC.",
                               "T.","NK.","Endothelial.","Fibroblast.","Epithelial.","Smooth_muscle.","Mast_3_gene."),search = T)

# Visualize gene signatures in violin plots: 
##########################################################################################################
# StackedVlnPlot(rna,features = c("B.1","Plasma.2","Mast.3","Macrophage.4",
#                                 "DC.5","T.6","NK.7","Endothelial.8","Fibroblast.9",
#                                 "Epithelial.10","Smooth_muscle.11"))+ggsave("Panglaodb_Signatures_Violin.pdf",width = 8,height = 16)
# 
# 
# StackedVlnPlot(rna,features = c("TPSB2","TPSAB1","KIT"))+ggsave("Mast_Signatures_Violin.pdf",width = 8,height = 4)

######################################################################################################
# Assess Mast cell enrichment to potentially rename clusters
# Assess Mast cell enrichment to potentially rename clusters
vln.df <- VlnPlot(rna,features = "Mast_3_gene.12")
library(psych)
data.mast <- describeBy(vln.df$data$Mast_3_gene.12, vln.df$data$ident, mat = TRUE)
data.mast <- dplyr::filter(data.mast,median > 0.225)

# Assess B cell enrichment to potentially rename clusters
vln.df <- VlnPlot(rna,features = "B.1")
library(psych)
data.B <- describeBy(vln.df$data$B.1, vln.df$data$ident, mat = TRUE)
data.B <- dplyr::filter(data.B,median > 0.225)



# Annotate mast/b cells
rna$mast.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.mast$group1),TRUE,FALSE)
rna$B.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.B$group1),TRUE,FALSE)


# Append SingleR annotations to cluster labels:
# The most common SingleR label in each cluster becomes the cluster label 

cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(SingleR) 


cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$SingleR[1]
  
}

# Rename cluster if median enrichment score is greater than 0.1  
if(nrow(data.mast) > 0){
  for (i in 1:nrow(data.mast)){
    cluster.ids[[data.mast$group1[i]]] <- "Mast cell"#Marker Mast cell cluster 
  }
}else{
  cluster.ids <- cluster.ids
}

# Rename cluster if median enrichment score is greater than 0.1  
if (nrow(data.B) >0 ){
  for (i in 1:nrow(data.B)){
    cluster.ids[[data.B$group1[i]]] <- "B cell"#Marker B cell cluster 
  }
}else{
  cluster.ids <- cluster.ids
}


cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$cell.type <- Idents(rna)

rna$cell.type <- paste0(rna$RNA_snn_res.0.7,"-",rna$cell.type)

Idents(rna) <- rna$cell.type
saveRDS(rna,"./endo_ovar_All_scRNA_processed-regressions.rds")

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cell.type <- rna$cell.type
#Manually annotate 24-cluster as smooth muscle
rna.df$cell.type <- str_replace_all(rna.df$cell.type,"24-Stromal fibroblast","24-Smooth muscle cells")

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

my_levels <- c(12,20,21,22,33,
               19,36,
               2,
               7,9,
               16,17,38,
               0,29,
               6,8,10,13,14,18,23,26,27,28,
               11,24,
               4,25,35,
               1,3,31,
               5,15,32,
               34,
               30,37)


# Relevel object@ident
rna.df$cluster.new <- factor(x = rna.df$cluster, levels = my_levels)


epithelial.cols <- colorRampPalette(c("#A0E989", "#245719"))
epithelial.cols <- epithelial.cols(15)

fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(10)

smooth.cols <- c("#b47fe5","#d8b7e8")

endo.cols <- c("#93CEFF","#4A99FF","#00478a")

t.cols <- c("gray60","gray20","gray80")

macro.cols <- c("#ff6600","#ff9d5c","#c24d00")

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
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
p1 <- LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)

rna$log10nCountRNA <- log10(rna$nCount_RNA)
p2 <- FeaturePlot(rna,features = "log10nCountRNA")
p3 <- FeaturePlot(rna,features = "nFeature_RNA")
CombinePlots(list(p1,p2,p3),ncol=3)+ggsave("Technical_Clusters_UMAP_RNA-regressions.pdf",width = 18,height = 6)
###########################################################################


# Plot nCountRNA and nFeatureRNA per tumor cluster
meta <- rna@meta.data
meta <- meta[meta$RNA_snn_res.0.7 %in% c(12,20,21,22,33,
                                         19,36,
                                         2,
                                         7,9,
                                         16,17,38,
                                         0,29),]

cells <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(meta$RNA_snn_res.0.7))))
names(cluster.ids) <- levels(factor(meta$RNA_snn_res.0.7))

for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)

levels(factor(meta$seurat_clusters))

cols <- test$Color

p1 <- ggplot(meta,aes(x=seurat_clusters,y=log10nCountRNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p2 <- ggplot(meta,aes(x=seurat_clusters,y=nFeature_RNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 500,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per epithelial cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scRNA_PercentMito_per_Cluster-epi.pdf",width = 8,height = 3)



# Plot nCountRNA and nFeatureRNA per fibroblast cluster
meta <- rna@meta.data
idx <- levels(factor(meta$cell.type))
idx <- grep("ibroblast",levels(factor(meta$cell.type)))
fibroblast <- levels(factor(meta$cell.type))[idx]
fibroblast <- fibroblast[-c(1,5,8,12,13)]# Remove non-normal fibroblasts
meta <- meta[meta$cell.type %in% fibroblast,]

cells <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(meta$RNA_snn_res.0.7))))
names(cluster.ids) <- levels(factor(meta$RNA_snn_res.0.7))

for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)

levels(factor(meta$seurat_clusters))

cols <- test$Color



p3 <- ggplot(meta,aes(x=seurat_clusters,y=log10nCountRNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p4 <- ggplot(meta,aes(x=seurat_clusters,y=nFeature_RNA))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 500,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per fibroblast cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scRNA_PercentMito_per_Cluster-fibro.pdf",width = 8,height = 3)

CombinePlots(list(p1,p2,p3,p4),ncol=1)+ggsave("BoxPlots_scRNA_variables-regressions.pdf",width = 6,height = 12)


# Patient Proportion per Cluster
# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c(12,20,21,22,33,
                                                      19,36,
                                                      2,
                                                      7,9,
                                                      16,17,38,
                                                      0,29,
                                                      6,8,10,13,14,18,23,26,27,28,
                                                      11,24,
                                                      4,25,35,
                                                      1,3,31,
                                                      5,15,32,
                                                      34,
                                                      30,37))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA-regressions.pdf",width = 4,height = 8)





###########################################################
# Part 3: ATAC clusters
###########################################################
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

# Make TSS UMAP
atac$barcode <- rownames(atac@cellColData)
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "cluster"
atac.df$cluster <- sub(".*?-","",atac.df$cluster)
atac.df$cluster <- sub("?-.*","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac$idx <- 1:nrow(atac@cellColData)

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

atac$cluster <- sub("?-.*","",atac$predictedGroup_ArchR)

length(which(atac$cluster==atac.df$cluster))
atac.df$TSSEnrichment <- as.numeric(atac$TSSEnrichment)

p2 <- ggplot(atac.df,aes(x = x,y=y,color = TSSEnrichment))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_gradient(low="lightgrey",high="blue")

# Make log10(nFrags) UMAP
atac$barcode <- rownames(atac@cellColData)
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "predicted.id"
atac.df$cluster <- sub(".*?-","",atac.df$predicted.id)
atac.df$cluster <- sub("?-.*","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac$idx <- 1:nrow(atac@cellColData)

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

atac$cluster <- sub("?-.*","",atac$predictedGroup_ArchR)

length(which(atac$cluster==atac.df$cluster))
atac.df$log10nFrags <- as.numeric(log10(atac$nFrags))
atac.df$Sample <- atac$Sample
atac.df$TSSEnrichment <- atac$TSSEnrichment

p3 <- ggplot(atac.df,aes(x = x,y=y,color = log10nFrags))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_gradient(low="lightgrey",high="blue")


CombinePlots(list(p1,p2,p3),ncol=3)+ggsave("Technical_Clusters_UMAP_ATAC.pdf",width = 18,height = 6)
###########################################################################


# Plot TSSEnrichment and log10nFrags per tumor cluster
atac.df.sub <- atac.df
atac.df.sub <- atac.df.sub[atac.df.sub$cluster %in% c(11,20,21,22,31,19,34,3,9,10,16,17,0,27),]
atac.df.sub$cluster <- factor(atac.df.sub$cluster,sort(as.numeric(levels(factor(atac.df.sub$cluster)))))
cells <- atac.df.sub %>% dplyr::group_by(cluster) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(atac.df.sub$cluster))))
names(cluster.ids) <- levels(factor(atac.df.sub$cluster))

for ( i in factor(cells$cluster)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(cluster ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)


cols <- test$Color

p1 <- ggplot(atac.df.sub,aes(x=cluster,y=log10nFrags))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p2 <- ggplot(atac.df.sub,aes(x=cluster,y=TSSEnrichment))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per tumor cluster:")+
  geom_hline(yintercept = 2,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per epithelial cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scatac_PercentMito_per_Cluster-epi.pdf",width = 8,height = 3)



# Plot TSS Enrichment and log10nFrags per fibroblast cluster
atac.df$cell.type <- sub(".*?-","",atac.df$predicted.id)
idx <- grep("ibroblast",levels(factor(atac.df$cell.type)))
fibroblast <- levels(factor(atac.df$cell.type))[idx]
fibroblast <- fibroblast[-c(1,5,7,11)]# Remove non-normal fibroblasts
atac.df <- atac.df[atac.df$cell.type %in% fibroblast,]
atac.df$cluster <- factor(atac.df$cluster,sort(as.numeric(levels(factor(atac.df$cluster)))))


cells <- atac.df %>% dplyr::group_by(cluster) %>% dplyr::count(Sample) 

cluster.ids <- rep("fill",length(levels(factor(atac.df$cluster))))
names(cluster.ids) <- levels(factor(atac.df$cluster))

for ( i in factor(cells$cluster)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(cluster ==i) %>% dplyr::arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$Sample[1]
  
}

cluster.ids <- as.data.frame(cluster.ids)
cluster.ids$cluster <- rownames(cluster.ids)
colnames(cluster.ids) <- c("Sample","cluster")
test <- merge(cluster.ids,sampleAnnot,by="Sample")
test$cluster <- as.numeric(test$cluster)
test <- dplyr::arrange(test,cluster)

levels(factor(atac.df$cluster))

cols <- test$Color



p3 <- ggplot(atac.df,aes(x=cluster,y=log10nFrags))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 3.0,linetype="dashed")+
  geom_hline(yintercept = 5.0,linetype="dashed")

p4 <- ggplot(atac.df,aes(x=cluster,y=TSSEnrichment))+geom_boxplot(fill=cols)+
  ggtitle("Technical variables per fibroblast cluster:")+
  geom_hline(yintercept = 2,linetype="dashed")

# ggplot(meta,aes(x=seurat_clusters,y=percent.mt))+geom_boxplot()+
#   ggtitle("Technical variables per fibroblast cluster:")+
#   geom_hline(yintercept = 20,linetype="dashed")+
#   ggsave("scatac_PercentMito_per_Cluster-fibro.pdf",width = 8,height = 3)

CombinePlots(list(p1,p2,p3,p4),ncol=1)+ggsave("BoxPlots_scATAC_variables.pdf",width = 6,height = 12)


# Patient Proportion per Cluster
# Patient proportion per subcluster in atac:

atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
colnames(atac.df)[3] <- "predicted.id"
atac.df$cluster <- sub(".*?-","",atac.df$predicted.id)
atac.df$cluster <- sub("?-.*","",atac.df$cluster)
atac.df$idx <- rownames(atac.df)
atac$idx <- 1:nrow(atac@cellColData)

atac.df <- dplyr::arrange(atac.df,as.numeric(idx))

atac$cluster <- sub("?-.*","",atac$predictedGroup_ArchR)

length(which(atac$cluster==atac.df$cluster))
atac.df$log10nFrags <- as.numeric(log10(atac$nFrags))
atac.df$Sample <- atac$Sample
atac.df$TSSEnrichment <- atac$TSSEnrichment


df <- atac.df %>% dplyr::group_by(cluster) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 

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
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height = 8)



# Save session info:
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
