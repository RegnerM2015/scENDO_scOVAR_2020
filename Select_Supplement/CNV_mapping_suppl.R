library(Seurat)
library(RColorBrewer)
library(ggplot2)

rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")

cancer.clusters <- c("0-Fibroblast","27-Fibroblast",
                     "3-Epithelial cell",
                     "31-Unciliated epithelia 1",
                     "22-Unciliated epithelia 2",
                     "21-Unciliated epithelia 1",
                     "11-Unciliated epithelia 1",
                     "19-Epithelial cell","34-Epithelial cell",
                     "20-Ciliated",
                     "9-Epithelial cell",
                     "10-Epithelial cell",
                     "17-Epithelial cell",
                     "16-Fibroblast")


# Color RNA UMAP according to main Figure 1:

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


cells <- colnames(rna[,rna$cell.type %in% cancer.clusters])

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
p1+ggsave("Full_Cohort_UMAP.pdf",width=8,height = 8)



# Patient 1 
rna.3533EL <- readRDS("./endo_3533EL_scRNA_processed.rds")
VlnPlot(rna.3533EL,features = "Total_CNVs")
VlnPlot(rna.3533EL,features="EPCAM")
VlnPlot(rna.3533EL,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3533EL[,rna.3533EL$RNA_snn_res.0.7 %in% c(10)]),"_1")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3533EL[,rna.3533EL$RNA_snn_res.0.7 %in% c(6)]),"_1")]
table(rna.sub$RNA_snn_res.0.7)

# Patient 2
rna.3571DL <- readRDS("./endo_3571DL_scRNA_processed.rds")
VlnPlot(rna.3571DL,features = "Total_CNVs")
VlnPlot(rna.3571DL,features="EPCAM")
VlnPlot(rna.3571DL,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3571DL[,rna.3571DL$RNA_snn_res.0.7 %in% c(4)]),"_2")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3571DL[,rna.3571DL$RNA_snn_res.0.7 %in% c(5)]),"_2")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3571DL[,rna.3571DL$RNA_snn_res.0.7 %in% c(9)]),"_2")]
table(rna.sub$RNA_snn_res.0.7)


# Patient 3
rna.36186L <- readRDS("./endo_36186L_scRNA_processed.rds")
VlnPlot(rna.36186L,features = "Total_CNVs")
VlnPlot(rna.36186L,features="EPCAM")
VlnPlot(rna.36186L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(0)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(1)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(2)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(3)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(4)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(7)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(8)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36186L[,rna.36186L$RNA_snn_res.0.7 %in% c(13)]),"_3")]
table(rna.sub$RNA_snn_res.0.7)

# Patient 4
rna.36639L <- readRDS("./endo_36639L_scRNA_processed.rds")
VlnPlot(rna.36639L,features = "Total_CNVs")
VlnPlot(rna.36639L,features="EPCAM")
VlnPlot(rna.36639L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36639L[,rna.36639L$RNA_snn_res.0.7 %in% c(1)]),"_4")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36639L[,rna.36639L$RNA_snn_res.0.7 %in% c(7)]),"_4")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.36639L[,rna.36639L$RNA_snn_res.0.7 %in% c(11)]),"_4")]
table(rna.sub$RNA_snn_res.0.7)

# Patient 5
rna.366C5L <- readRDS("./endo_366C5L_scRNA_processed.rds")
VlnPlot(rna.366C5L,features = "Total_CNVs")
VlnPlot(rna.366C5L,features="EPCAM")
VlnPlot(rna.366C5L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.366C5L[,rna.366C5L$RNA_snn_res.0.7 %in% c(10)]),"_5")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.366C5L[,rna.366C5L$RNA_snn_res.0.7 %in% c(16)]),"_5")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.366C5L[,rna.366C5L$RNA_snn_res.0.7 %in% c(17)]),"_5")]
table(rna.sub$RNA_snn_res.0.7)


# Patient 6
rna.37EACL <- readRDS("./endo_ovar_metastasis_37EACL-ovarian_reference_scRNA_processed.rds")
VlnPlot(rna.37EACL,features = "Total_CNVs")
VlnPlot(rna.37EACL,features="EPCAM")
VlnPlot(rna.37EACL,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.37EACL[,rna.37EACL$RNA_snn_res.0.7 %in% c(0)]),"_6")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.37EACL[,rna.37EACL$RNA_snn_res.0.7 %in% c(9)]),"_6")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.37EACL[,rna.37EACL$RNA_snn_res.0.7 %in% c(13)]),"_6")]
table(rna.sub$RNA_snn_res.0.7)



# Patient 7
rna.38FE7L <- readRDS("./ovar_38FE7L_scRNA_processed.rds")
VlnPlot(rna.38FE7L,features = "Total_CNVs")
VlnPlot(rna.38FE7L,features="EPCAM")
VlnPlot(rna.38FE7L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(0)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(1)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(3)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(4)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(7)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(9)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(10)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.38FE7L[,rna.38FE7L$RNA_snn_res.0.7 %in% c(11)]),"_7")]
table(rna.sub$RNA_snn_res.0.7)



# Patient 8
rna.3BAE2L <- readRDS("./ovar_3BAE2L_scRNA_processed.rds")
VlnPlot(rna.3BAE2L,features = "Total_CNVs")
VlnPlot(rna.3BAE2L,features="EPCAM")
VlnPlot(rna.3BAE2L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(0)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(1)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(7)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(9)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(11)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(14)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3BAE2L[,rna.3BAE2L$RNA_snn_res.0.7 %in% c(18)]),"_8")]
table(rna.sub$RNA_snn_res.0.7)



# Patient 9
rna.3E5CFL <- readRDS("./ovar_3E5CFL_scRNA_processed.rds")
VlnPlot(rna.3E5CFL,features = "Total_CNVs")
VlnPlot(rna.3E5CFL,features="EPCAM")
VlnPlot(rna.3E5CFL,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E5CFL[,rna.3E5CFL$RNA_snn_res.0.7 %in% c(0)]),"_11")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E5CFL[,rna.3E5CFL$RNA_snn_res.0.7 %in% c(2)]),"_11")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E5CFL[,rna.3E5CFL$RNA_snn_res.0.7 %in% c(5)]),"_11")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E5CFL[,rna.3E5CFL$RNA_snn_res.0.7 %in% c(6)]),"_11")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E5CFL[,rna.3E5CFL$RNA_snn_res.0.7 %in% c(9)]),"_11")]
table(rna.sub$RNA_snn_res.0.7)



# Patient 10
rna.3CCF1L <- readRDS("./ovar_3CCF1L_scRNA_processed.rds")
VlnPlot(rna.3CCF1L,features = "Total_CNVs")
VlnPlot(rna.3CCF1L,features="EPCAM")
VlnPlot(rna.3CCF1L,features="MUC16")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(0)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(1)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(7)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(8)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(10)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(13)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(14)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3CCF1L[,rna.3CCF1L$RNA_snn_res.0.7 %in% c(15)]),"_9")]
table(rna.sub$RNA_snn_res.0.7)




# Patient 11
rna.3E4D1L <- readRDS("./ovar_3E4D1L_scRNA_processed.rds")
VlnPlot(rna.3E4D1L,features = "Total_CNVs")
VlnPlot(rna.3E4D1L,features = "TOP2A")
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(0)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(2)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(4)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(5)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(8)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(10)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
rna.sub <- rna[,colnames(rna) %in% paste0(colnames(rna.3E4D1L[,rna.3E4D1L$RNA_snn_res.0.7 %in% c(1)]),"_10")]
table(rna.sub$RNA_snn_res.0.7)
