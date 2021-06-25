###################################################
# Matt Regner
# Franco Lab 
# Nov 2020
# Description:Perform GSVA with CancerSEA signatures
###################################################


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




source("./stacked_violin.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(stringr)
library(ComplexHeatmap)


# Read in gene sets and store in gmt format:
files <- list.files(pattern = ".txt")

names <- str_remove(files,".txt")

desired_length <- 14 # or whatever length you want
filler <- vector(mode = "list", length = desired_length)
names(filler) <- names
for (i in names(filler)){
  genes <- read.delim(paste0(i,".txt"),header = T)
  
  genes.vector <- genes$GeneName
  
  filler[[i]] <- genes.vector
  
}

gset <- filler
# Set labels of cell types of interest




# Using RNA clusters that are likely malignant 
labels <- c("3-Epithelial cell",
            "9-Epithelial cell",
            "10-Epithelial cell",
            "11-Unciliated epithelia 1",
            "17-Epithelial cell",
            "19-Epithelial cell" ,
            "21-Unciliated epithelia 1",
            "22-Unciliated epithelia 2",
            "31-Unciliated epithelia 1",
            "34-Epithelial cell",
            "20-Ciliated" ,
            "16-Fibroblast",
            "0-Fibroblast" ,
            "27-Fibroblast")# Screen for "pure" patient clusters

# Pseudobulk RNA expression
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")

FeaturePlot(rna,features = "Total_CNVs")+ggsave("CNV.pdf")

rna.counts <- rna@assays$RNA@counts

rna.pseudobulk <- data.frame(rownames(rna))
for (i in labels){
  cells <- rownames(dplyr::filter(rna@meta.data,cell.type == i))
  
  rna.counts.sub <- rna.counts[,colnames(rna.counts) %in% cells]
  
  rna.counts.bulk <- rowSums(as.matrix(rna.counts.sub))
  
  rna.pseudobulk$i <- rna.counts.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
}
rownames(rna.pseudobulk) <- rna.pseudobulk[,1]
rna.pseudobulk <- rna.pseudobulk[,-1]
dim(rna.pseudobulk)


res <- GSVA::gsva(as.matrix(rna.pseudobulk),gset.idx.list = gset,method = "gsva",kcdf = "Poisson")
head(res)

ha = HeatmapAnnotation(
  Site = factor(c(rep("1",7),rep("2",7))),
  Type = factor(c(rep("1",7),rep("2",3),rep("3",4))),
  Histology = factor(c(rep("1",7),rep("2",3),rep("3",2),rep("4",2))),
  Stage = factor(c(rep("1",4),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:14)))
)
pdf("GSVA_Heatmap_raw_poisson.pdf",width = 6,height = 8)
# Make heatmap annotation
Heatmap(scale(res),cluster_rows = T,top_annotation = ha)
dev.off()



# Rlog transform raw counts
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.rna <- rlog(dds,blind = T)
dat <- as.matrix(assay(rlog.rna))


res <- GSVA::gsva(dat,gset.idx.list = gset,method = "gsva",kcdf = "Gaussian")
head(res)

                  
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",7),rep("2",7))),
  Type = factor(c(rep("1",7),rep("2",3),rep("3",4))),
  Histology = factor(c(rep("1",7),rep("2",3),rep("3",2),rep("4",2))),
  Stage = factor(c(rep("1",4),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:14)))
)
pdf("GSVA_Heatmap_Rlog_gaussian.pdf",width = 6,height =10)
# Make heatmap annotation
Heatmap(scale(res),cluster_rows = T,top_annotation = ha)
dev.off()



# Make UMAP plot highlighting by cell clusters of interest:
rna$gsva.study <- ifelse(rna$cell.type %in% labels,TRUE,FALSE)
p1 <-DimPlot(rna,group.by = "gsva.study",cols = c("gray48","forestgreen"))+ggsave("GSVA_Study_UMAP.pdf",width = 6,height = 4)


# Stacked violin for all GSVA signatures in AddModuleScore()

rna.sub <- subset(rna, idents = labels)
rna.sub <- AddModuleScore(rna.sub,features = gset,search = T,names = names)

Idents(rna.sub) <- "RNA_snn_res.0.7"

my_levels <- rev(c(11,20,21,22,31,19,34,
                   3,9,10,17,16,
                   0,27))

# Relevel object@ident
rna.sub$cluster.new <- factor(x = rna.sub$RNA_snn_res.0.7, levels = my_levels)
Idents(rna.sub) <- "cluster.new"


StackedVlnPlot(rna.sub,features = c("Cluster9",
                                "Cluster10",
                                "Cluster7"))+ggsave("CancerSEA_Vln.pdf",width = 12,height =10)

# Check significance for each Gene signature (3,7,10)
summary(aov(Cluster7 ~ cluster.new, data = rna.sub@meta.data))
kruskal.test(Cluster7 ~ cluster.new, data = rna.sub@meta.data)

summary(aov(Cluster10 ~ cluster.new, data = rna.sub@meta.data))
kruskal.test(Cluster10 ~ cluster.new, data = rna.sub@meta.data)

summary(aov(Cluster9 ~ cluster.new, data = rna.sub@meta.data))
kruskal.test(Cluster9 ~ cluster.new, data = rna.sub@meta.data)

# Show proportion of patient cells per cluster:
meta <- rna.sub@meta.data

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 

df %>% 
  dplyr::mutate(cluster.new = factor(Cluster,levels = c(11,20,21,22,31,19,34,
                                                        3,9,10,17,16,
                                                        0,27))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cluster.new))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 6)




writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

