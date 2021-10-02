###################################################
# Matt Regner
# Franco Lab 
###################################################

library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(ArchR)
library(stringr)
library(stringi)
library(RColorBrewer)
library(dplyr)
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



# Read in RNA data for full cohort
rna <- readRDS("endo_ovar_All_scRNA_processed.rds")


# Read in matching ATAC data for full cohort (ArchR project after label transfer)
atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")




# Plot Patient UMAPs for RNA/ATAC

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$Sample <- rna$Sample

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Sample_RNA.pdf",width = 8,height = 7)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$Sample <- gsub(".*-","",atac.df$color)

ggplot(atac.df,aes_string(x = "x",y="y",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 7)




# Plot cell type UMAPs for RNA/ATAC

# RNA
levels(factor(rna$RNA_snn_res.0.7))
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


# ATAC
levels(factor(atac$predictedGroup_ArchR))

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

p1 <- ggplot(atac.df,aes(x =x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_ATAC.pdf",width = 8,height = 7)
p1 <- ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_ATAC-labels.pdf",width = 8,height = 7)

#########################################################################################################

# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
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


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% group_by(predictedGroup_ArchR) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
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



# Write cluster total # of cells to output files
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")


total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")


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

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
