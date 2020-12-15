###################################################
# Matt Regner
# Franco Lab 
# Nov 2020
# Description: plot patient UMAPs, cell type UMAPs,
# histology UMAPs, patient proportion per cluster
###################################################
source("./P2G_Heatmap_Distal.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(ArchR)
library(stringr)
library(stringi)
library(dplyr)
library(bedtoolsr)
library(tidyr)
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
rna.sample.plot +ggsave("Sample_RNA.pdf",width = 8,height = 6)

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
ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)




# Plot Histology UMAPs for RNA/ATAC

# Add Histology variable to RNA data
rna.df$Histology <- rna.df$Sample
rna.df$Histology <- str_replace_all(rna.df$Histology, "3533EL", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3571DL", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "36186L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "36639L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "366C5L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "37EACL", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "38FE7L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3BAE2L", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3E5CFL", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3E4D1L", "GIST")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3CCF1L", "Carcinosarcoma")


rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Histology))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = c("gray60","coral2","mediumpurple1","#00CED1"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Histology_RNA.pdf",width = 8,height = 6)


# Add Histology variable to ATAC data
atac.df$Histology <- atac.df$Sample
atac.df$Histology <- str_replace_all(atac.df$Histology, "3533EL", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3571DL", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "36186L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "36639L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "366C5L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "37EACL", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "38FE7L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3BAE2L", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3E5CFL", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3E4D1L", "GIST")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3CCF1L", "Carcinosarcoma")


atac.sample.plot <-ggplot(atac.df,aes(x = x,y=y,color = Histology))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = c("gray60","coral2","mediumpurple1","#00CED1"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
atac.sample.plot +ggsave("Histology_ATAC.pdf",width = 8,height = 6)





# Plot cell type UMAPs for RNA/ATAC

# RNA
levels(factor(rna$cell.type))
rna.df$cell.type <- rna$RNA_snn_res.0.7


my_levels <-c(1,9,10,11,18,19,21,22,23,32,35,0,7,8,12,14,15,16,17,24,25,26,27,29,30,4,6,20,34,2,3,31,5,13,33,28,36
)

# Relevel object@ident
rna.df$cluster.new <- factor(x = rna.df$cell.type, levels = my_levels)


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(11)

fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(14)

smooth.cols <- c("#b47fe5")

endo.cols <- c("#93CEFF","#4A99FF","#01579B")

t.cols <- c("gray60","gray20","gray40")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)

ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_RNA.pdf",width = 8,height = 6)



# ATAC
levels(factor(atac$predictedGroup_ArchR))


atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- atac.df$data

atac.df$cell.type <- sub(".*?-","",atac.df$color)
atac.df$cell.type <- as.factor(atac.df$cell.type)
atac.df$cell.type <- gsub("-.*","",atac.df$cell.type)


my_levels <-c(1,9,10,11,18,19,21,22,23,32,35,0,7,8,12,14,15,16,17,24,25,26,27,29,30,4,6,20,34,2,3,31,5,13,33,28,36
)

# Relevel object@ident
atac.df$cluster.new <- factor(x = atac.df$cell.type, levels = my_levels)


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(11)

fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(14)

smooth.cols <- c("#b47fe5")

endo.cols <- c("#93CEFF","#4A99FF","#01579B")

t.cols <- c("gray60","gray20","gray40")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)

ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_ATAC.pdf",width = 8,height = 6)


#########################################################################################################

# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("1","9","10","11","18","19","21","22","23","32","35",
                                                      "0","7","8","12","14","15","16","17","24","25","26","27","29","30",
                                                      "4",
                                                      "6","20","34",
                                                      "2","3","31",
                                                      "5","13",
                                                      "33",
                                                      "28","36"))) %>% 
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
  dplyr::mutate(cell.type = factor(Cluster,levels = c("1","9","10","11","18","19","21","22","23","32","35",
                                                      "0","7","8","12","14","15","16","17","24","25","26","27","29","30",
                                                      "4",
                                                      "6","20","34",
                                                      "2","3","31",
                                                      "5","13",
                                                      "33",
                                                      "28","36"))) %>% 
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





# Plot peak2gene heatmap
###########################################################
atac.peaks <- readRDS("./final_archr_proj_archrGS.rds")

sig.distal <- readRDS("./Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")
sig.distal <- dplyr::filter(sig.distal,EmpFDR <= 0.05)

full <- data.frame(chrom = sig.distal$seqnames.x,start = sig.distal$start.x,end = sig.distal$end.x)
full.gr <- makeGRangesFromDataFrame(full)
empirical.sig.All.peaks.from.full.cohort <- makeGRangesFromDataFrame(full.gr)
###################################################################################

# Subset colors to new cell cluster labels

atac.peaks@cellColData$cluster.num <- gsub("-.*","",atac.peaks@cellColData$predictedGroup_ArchR)

my_levels <-c(1,9,10,11,18,19,21,22,23,32,35,0,7,8,12,14,15,16,17,24,25,26,27,29,30,4,6,20,34,2,3,31,5,13,33,28,36
  )

# Relevel object@ident
atac.peaks@cellColData$cluster.new <- factor(x = as.numeric(atac.peaks@cellColData$cluster.num), levels = my_levels)
#atac.peaks@cellColData$cluster.new <- as.character(atac.peaks@cellColData$cluster.new)
names(cols) <- levels(factor(atac.peaks@cellColData$cluster.new))


# For every gene in user defined list, plot P2G linkage browser track with 
# 1) All ENCODE cCREs
# 2) Distal ENCODE cCREs
# 3) All statistically significant peaks that overlap with EmpFDR significant peaks from full cohort
# 4) All statistically signiicant peaks (EmpFDR is calculated within All cohort only)
# Next, plot matching scRNA-seq expression in violin plot
########################################################################################################

# Read in ENCODE cCRE information:
encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
encode.all <- makeGRangesFromDataFrame(encode.all)

encode.distal <- read.delim("./GRCh38-ccREs.dELS.bed",header =F)
colnames(encode.distal)[1:3] <- c("seqnames","start","end")
encode.distal <- makeGRangesFromDataFrame(encode.distal)

# Read in ENCODE DN-ase epithelium 
encode.epithelium <- read.delim("./ENCODE_Epithelium_DNase.bed",header = T,sep = "\t")
colnames(encode.epithelium)[1:3] <- c("seqnames","start","end")
encode.epithelium <- makeGRangesFromDataFrame(encode.epithelium)





# Annotation 1:
sig.peaks <- readRDS("./P2G_Hits_RHEB_data.rds")
sig.peaks <- as.data.frame(sig.peaks$peak)
sig.peaks$start <- sig.peaks$`sig.peaks$peak`
sig.peaks$end <- sig.peaks$`sig.peaks$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)


       

levels(factor(atac.peaks$predictedGroup_ArchR))
levels(factor(atac.peaks$cluster.new))
levels(factor(atac.peaks@cellColData$cluster.new))# Order that we want 
cols# also order that we want 

# library(plyr)
# atac.peaks@cellColData$cluster.new <- revalue(atac.peaks@cellColData$cluster.new, c(
#   "1"="1-1", 
#   "9"="2-9",
#   "10"="3-10",
#   "11"="4-11",
#   "18"="5-18",
#   "19"="6-19",
#   "21"="7-21",
#   "22"="8-22",
#   "23"="9-23",
#   "32"="10-32",
#   "35"="11-35",
#   "0"="12-0",
#   "7"="13-7",
#   "8"="14-8",
#   "12"="15-12",
#   "14"="16-14",
#   "15"="17-15",
#   "16"="18-16",
#   "17"="19-17",
#   "24"="20-24",
#   "25"="21-25",
#   "26"="22-26",
#   "27"="23-27",
#   "29"="24-29",
#   "30"="25-30",
#   "4"="26-4",
#   "6"="27-6",
#   "20"="28-20",
#   "34"="29-34",
#   "2"="30-2",
#   "3"="31-3",
#   "31"="32-31",
#   "5"="33-5",
#   "13"="34-13",
#   "33"="35-33",
#   "28"="36-28",
#   "36"="37-36"
#   ))

plot <- plotBrowserTrack(atac.peaks,geneSymbol ="RHEB", groupBy = "cluster.new",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium), 
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001),pal=cols,upstream = 20000,downstream = 50000)

pdf("RHEB_final.pdf",width = 6,height = 12)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="PRKAG2", groupBy = "cluster.new",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.All.peaks.from.full.cohort,TrackE=sig.peaks), 
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001),pal=cols,upstream = 500000,downstream = 20000)

pdf("RHEB_final_sig_peaks_anno.pdf",width = 6,height = 12)
grid::grid.draw(plot[[1]])
dev.off()




# Differential peak check:

markersPeaks <- getMarkerFeatures(
  ArchRProj = atac.peaks,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList.atac <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerList.atac <- as.data.frame(markerList.atac)
markerList.atac$peak <- paste0(markerList.atac$seqnames,":",markerList.atac$start,"-",markerList.atac$end)




# Make violin plots for RHEB expression and mTOR pathway expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "RHEB",pt.size = 0)+coord_flip()+NoLegend()+ggsave("RHEB_vln_new.pdf",width = 4,height = 8)

p1 <- p1$data
colnames(p1)
kruskal.test(RHEB~ident,data = p1)


gset <- read.delim("./BIOCARTA_MTOR_PATHWAY.txt",header = T)
gset <- as.character(gset$BIOCARTA_MTOR_PATHWAY[2:length(gset$BIOCARTA_MTOR_PATHWAY)])
rna <- AddModuleScore(rna,features = list(gset),name = "mTOR_members",search = T)


p1 <- VlnPlot(rna,features = "mTOR_members1",pt.size = 0)+coord_flip()+NoLegend()+
  ggsave("mTOR_pathway_vln_new.pdf",width = 4,height = 8)

p1 <- p1$data
colnames(p1)
kruskal.test(mTOR_members1~ident,data = p1)
