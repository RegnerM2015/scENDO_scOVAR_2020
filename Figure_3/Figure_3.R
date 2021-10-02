###########################################################
# Matt Regner
# Franco Lab
# Date: 2020-2021
#
# Description: Figure 3, EEC enhancer discovery 
# RNA and ATAC data
###########################################################
source("./P2G_Heatmap_Distal.R")
source("./Archr_Peak_Null_Permute.R")
source("./ArchRBrowser.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(bedtoolsr)
library(ArchR)
library(stringr)
library(stringi)

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


###########################################################
# Read in RNA data for full cohort
rna <- readRDS("endo_EEC_scRNA_processed.rds")


# Read in matching ATAC data for full cohort (ArchR project after label transfer)
atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")
###########################################################

# Plot Patient UMAPs for RNA/ATAC
###########################################################
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
  scale_color_manual(values = sampleColors[1:5])+
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
  scale_color_manual(values = sampleColors[1:5])+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)
###########################################################




# Plot cell type UMAPs for RNA/ATAC
###########################################################
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
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_RNA-labels.pdf",width = 8,height = 7)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)


atac.df$cell.type <- str_replace(atac.df$cell.type,pattern = "12-Stromal fibroblasts","12-Smooth muscle cells")
#Rename 12 cluster as smooth muscle based on Myh11 expression
atac.df$cluster <- sub("?-.*","",atac.df$cell.type)
atac.df$cluster.new  <- factor(atac.df$cluster,levels = c("6",
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

ggplot(atac.df,aes_string(x = "x",y="y",color = "cluster.new"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("Cell_type_ATAC.pdf"),width = 8,height = 6)

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

###########################################################

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
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_RNA.pdf",width =4,height = 11)


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% group_by(predictedGroup_ArchR) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
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
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height =11)
###########################################################


# Write cluster total # of cells to output files
###########################################################
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")


total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")
###########################################################


# Suppl. Figure CNV plot
# Plot CNV box: 
meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7,levels = rev(c("6",
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
                                                           "27","28")))

ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot()+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  ggsave("CNV_BoxPlot.pdf",width = 4,height = 8)

# clear up memory:
rm(atac)
rm(rna)


# Make modified getP2G function:
getPeak2GeneLinks.mod <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  PValCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = PValCutOff, name = "PValCutOff", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    return(NULL)
  }
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    
    return(NULL)
    
  }else{
    
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$RawPVal <= PValCutOff), ,drop=FALSE]
    
    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    
    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }
    
    if(returnLoops){
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(Peak2GeneLinks = loops)
      
      return(loops)
      
    }else{
      
      return(p2g)
      
    }
    
  }
  
}

# Read in other annotation features:
cancer.specific <- readRDS("Cancer_specific_P2G_table.rds")



atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
# ATAC
levels(factor(atac$predictedGroup_ArchR))
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
               "23",
               "4" ,"5","12",
               "3","11","26",
               "0","18",
               "24","7",
               "25",
               "27","28")

for ( i in levels(factor(atac$predictedGroup_ArchR))){
  num <-  gsub("-.*","",i)
  idx <- match(num,my_levels)
  atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = i,replacement = paste0(idx,"_",atac$predictedGroup_ArchR))
  print("iter complete")
}


encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
encode.all <- makeGRangesFromDataFrame(encode.all)

ft.peaks <- readRDS("./Fallopian_Tube_Cell_line_Peaks.rds")
ft.peaks <- ft.peaks[,3:5]
colnames(ft.peaks)[1:3] <- c("seqnames","start","end")
ft.peaks <- makeGRangesFromDataFrame(ft.peaks)

ov.peaks <- readRDS("./Ovarian_Epithelial_Cell_line_Peaks.rds")
ov.peaks <- ov.peaks[,3:5]
colnames(ov.peaks)[1:3] <- c("seqnames","start","end")
ov.peaks <- makeGRangesFromDataFrame(ov.peaks)
# 28 B cell has too few cells to plot
plot <- plotBrowserTrack(atac,geneSymbol = "SOX9", groupBy = "predictedGroup_ArchR",
                         loops = getPeak2GeneLinks(atac,FDRCutOff = 1,varCutOffATAC = 0,
                                                   varCutOffRNA = 0),upstream = 9000,downstream =202000,
                         features=GRangesList(TrackA=encode.all,TrackB=ft.peaks,TrackC=ov.peaks),
                         pal= cols[-c(15,29)],
                         ylim=.9965)

pdf("SOX9_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

plot <- plotBrowserTrack(atac,geneSymbol = "CD24", groupBy = "predictedGroup_ArchR",
                         loops = getPeak2GeneLinks(atac,FDRCutOff = 1,varCutOffATAC = 0,
                                                   varCutOffRNA = 0),upstream =110000,downstream =3000,
                         features=GRangesList(TrackA=encode.all,TrackB=ft.peaks,TrackC=ov.peaks),
                         pal= cols[-c(15,29)],
                         ylim=.996)

pdf("CD24_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

plot <- plotBrowserTrack(atac,geneSymbol = "IMPA2", groupBy = "predictedGroup_ArchR",
                         loops = getPeak2GeneLinks(atac,FDRCutOff = 1,varCutOffATAC = 0,
                                                   varCutOffRNA = 0),upstream =68000,downstream =12500,
                         features=GRangesList(TrackA=encode.all,TrackB=ft.peaks,TrackC=ov.peaks),
                         pal= cols[-c(15,29)],
                         ylim=.9910)

pdf("IMPA2_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

# Note: excluding 28-B cell b/c it has too few cells in ATAC browser track
# Plot matching scRNA-seq data:
rna <- readRDS("endo_EEC_scRNA_processed.rds")

my_levels <- rev(c("6",
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
                   "23",
                   "4" ,"5","12",
                   "3","11","26",
                   "0","18",
                   "24","7",
                   "25",
                   "27"))


# Relevel object@ident
rna.sub <- rna[,rna$RNA_snn_res.0.7 %in% my_levels]
rna.sub@active.ident <- factor(x =rna.sub$RNA_snn_res.0.7, levels =my_levels)


p1 <- VlnPlot(rna.sub,features = "SOX9",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=SOX9))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("SOX9 expression")+ggsave("Vln.pdf",width=3,height = 8)


# Check SOX9 statistical significance:
test <- FindMarkers(rna.sub,ident.1 = c("6",
                                        "8",
                                        "10",
                                        "14",
                                        "15",
                                        "19",
                                        "21",
                                        "22"),
                    ident.2 = c("1",
                                "2",
                                "9",
                                "13",
                                "16",
                                "17",
                                "23",
                                "4" ,"5","12",
                                "3","11","26",
                                "0","18",
                                "24","7",
                                "25",
                                "27"),only.pos = T)
test$genes <- rownames(test)
test <- test[test$gene == "SOX9",]
print(test)

p1 <- VlnPlot(rna.sub,features = "CD24",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=CD24))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CD24 expression")+ggsave("Vln.pdf",width=3,height = 8)


# Check CD24 statistical significance:
test <- FindMarkers(rna.sub,ident.1 = c("6",
                                        "8",
                                        "10",
                                        "14",
                                        "15",
                                        "19",
                                        "21",
                                        "22"),
                    ident.2 = c("1",
                                "2",
                                "9",
                                "13",
                                "16",
                                "17",
                                "23",
                                "4" ,"5","12",
                                "3","11","26",
                                "0","18",
                                "24","7",
                                "25",
                                "27"),only.pos = T)
test$genes <- rownames(test)
test <- test[test$gene == "CD24",]
print(test)


p1 <- VlnPlot(rna.sub,features = "IMPA2",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=IMPA2))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("IMPA2 expression")+ggsave("Vln-IMPA2.pdf",width=3,height = 8)


# Check IMPA2 statistical significance:
test <- FindMarkers(rna.sub,ident.1 = c("6",
                                        "8",
                                        "10",
                                        "14",
                                        "15",
                                        "19",
                                        "21",
                                        "22"),
                    ident.2 = c("1",
                                "2",
                                "9",
                                "13",
                                "16",
                                "17",
                                "23",
                                "4" ,"5","12",
                                "3","11","26",
                                "0","18",
                                "24","7",
                                "25",
                                "27"),only.pos = T,features = "IMPA2",logfc.threshold = 0)
print(test)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

