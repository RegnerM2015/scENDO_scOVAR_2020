library(ggplot2)
library(tidyverse)
library(ArchR)
library(dplyr)

p2g <- readRDS("./FullCohort_All_P2G_Observed.rds")
p2g <- p2g[p2g$Correlation >= 0.45,]
p2g <- p2g[p2g$RawPVal <= 1e-12,]
p2g <- p2g[p2g$peakType == "Distal",]# Subset to distal P2Gs
p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
p2g.cancer <- readRDS("./FullCohort_Cancer_specific_P2G_table.rds")
p2g$Cancer.Specific <- ifelse(p2g$idx %in% p2g.cancer$idx,"Cancer","Normal")


pdf("./Genes_Peaks_Histrograms-FullCohort.pdf",width=6,height =6)
par(mfrow=c(2,2)) 
hist(table(p2g$peakName),breaks=28)
abline(v=mean(table(p2g$peakName)),col="red")
hist(table(p2g$geneName),breaks=28)
abline(v=mean(table(p2g$geneName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName)),col="red")
dev.off()
# Number of peaks per number of genes
p2g.cancer <- readRDS("./FullCohort_Cancer_specific_P2G_table.rds")
p2g.normal <- readRDS("./FullCohort_All_P2G_Observed.rds")
p2g.normal <- p2g.normal[p2g.normal$Correlation >= 0.45,]
p2g.normal <- p2g.normal[p2g.normal$RawPVal <= 1e-12,]
p2g.normal <- p2g.normal[p2g.normal$peakType == "Distal",]
p2g.normal <- p2g.normal[p2g.normal$peakName %ni% p2g.cancer$peakName,]

df.cancer <- data.frame(num.genes = table(p2g.cancer$peakName))
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq < 3,"1-2",df.cancer$num.genes.Freq)
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq >2,"3+",df.cancer$cat)

df.normal <- data.frame(num.genes = table(p2g.normal$peakName))
df.normal$cat <- ifelse(df.normal$num.genes.Freq < 3,"1-2",df.normal$num.genes.Freq)
df.normal$cat <- ifelse(df.normal$num.genes.Freq >2,"3+",df.normal$cat)

head(df.normal)
head(df.cancer)

df.normal$type <- "normal"
df.cancer$type <- "cancer"

df<- rbind(df.normal,df.cancer)


p1 <- ggplot(df, aes(x=type, y=num.genes.Freq, fill = type)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(1.44,1.65))+NoLegend()

df.cancer %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.cancer)*100)) -> cancer
df.normal %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.normal)*100)) -> normal

cancer$type <- "cancer"
normal$type <- "normal"

comb <- rbind(cancer,normal)
print(comb)
p2 <- ggplot(comb,aes(x=cat,y=perc,fill=type))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=round(perc,3)), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+NoLegend()


CombinePlots(list(p2,p1),ncol=2)+ggsave("Barcharts_P2G_plots-FullCohort.pdf",width=6,height=3)



# Cancer specific peaks link to more genes on average with statistical significance:
test <- wilcox.test(df.cancer$num.genes.Freq,df.normal$num.genes.Freq,correct = T)
print(test)
print(test$p.value)






p2g <- readRDS("./EEC_All_P2G_Observed.rds")
p2g <- p2g[p2g$Correlation >= 0.45,]
p2g <- p2g[p2g$RawPVal <= 1e-12,]
p2g <- p2g[p2g$peakType == "Distal",]# Subset to distal P2Gs
p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
p2g.cancer <- readRDS("./EEC_Cancer_specific_P2G_table.rds")
p2g$Cancer.Specific <- ifelse(p2g$idx %in% p2g.cancer$idx,"Cancer","Normal")


pdf("./Genes_Peaks_Histrograms-EEC.pdf",width=6,height =6)
par(mfrow=c(2,2)) 
hist(table(p2g$peakName))
abline(v=mean(table(p2g$peakName)),col="red")
hist(table(p2g$geneName))
abline(v=mean(table(p2g$geneName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName)),col="red")
dev.off()
# Number of peaks per number of genes
p2g.cancer <- readRDS("./EEC_Cancer_specific_P2G_table.rds")
p2g.normal <- readRDS("./EEC_All_P2G_Observed.rds")

p2g.normal <- p2g.normal[p2g.normal$RawPVal <= 1e-12,]
p2g.normal <- p2g.normal[p2g.normal$Correlation >= 0.45,]
p2g.normal <- p2g.normal[p2g.normal$peakType == "Distal",]
p2g.normal <- p2g.normal[p2g.normal$peakName %ni% p2g.cancer$peakName,]

df.cancer <- data.frame(num.genes = table(p2g.cancer$peakName))
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq < 3,"1-2",df.cancer$num.genes.Freq)
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq >2,"3+",df.cancer$cat)

df.normal <- data.frame(num.genes = table(p2g.normal$peakName))
df.normal$cat <- ifelse(df.normal$num.genes.Freq < 3,"1-2",df.normal$num.genes.Freq)
df.normal$cat <- ifelse(df.normal$num.genes.Freq >2,"3+",df.normal$cat)

head(df.normal)
head(df.cancer)

df.normal$type <- "normal"
df.cancer$type <- "cancer"

df<- rbind(df.normal,df.cancer)


p1 <- ggplot(df, aes(x=type, y=num.genes.Freq, fill = type)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(1.35,1.62))+NoLegend()

df.cancer %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.cancer)*100)) -> cancer
df.normal %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.normal)*100)) -> normal

cancer$type <- "cancer"
normal$type <- "normal"

comb <- rbind(cancer,normal)
print(comb)
p2 <- ggplot(comb,aes(x=cat,y=perc,fill=type))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=round(perc,3)), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+NoLegend()


CombinePlots(list(p2,p1),ncol=2)+ggsave("Barcharts_P2G_plots-EEC.pdf",width=6,height=3)



# Cancer specific peaks link to more genes on average with statistical significance:
test <- wilcox.test(df.cancer$num.genes.Freq,df.normal$num.genes.Freq,correct = T)
print(test)
print(test$p.value)






p2g <- readRDS("./HGSOC_All_P2G_Observed.rds")
p2g <- p2g[p2g$RawPVal <= 1e-12,]
p2g <- p2g[p2g$Correlation >= 0.45,]

p2g <- p2g[p2g$peakType == "Distal",]# Subset to distal P2Gs
p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
p2g.cancer <- readRDS("./HGSOC_Cancer_specific_P2G_table.rds")
p2g$Cancer.Specific <- ifelse(p2g$idx %in% p2g.cancer$idx,"Cancer","Normal")


pdf("./Genes_Peaks_Histrograms-HGSOC.pdf",width=6,height =6)
par(mfrow=c(2,2)) 
hist(table(p2g$peakName))
abline(v=mean(table(p2g$peakName)),col="red")
hist(table(p2g$geneName))
abline(v=mean(table(p2g$geneName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Cancer",]$peakName)),col="red")
hist(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName))
abline(v=mean(table(p2g[p2g$Cancer.Specific == "Normal",]$peakName)),col="red")
dev.off()
# Number of peaks per number of genes
p2g.cancer <- readRDS("./HGSOC_Cancer_specific_P2G_table.rds")
p2g.normal <- readRDS("./HGSOC_All_P2G_Observed.rds")

p2g.normal <- p2g.normal[p2g.normal$RawPVal <= 1e-12,]
p2g.normal <- p2g.normal[p2g.normal$Correlation >= 0.45,]
p2g.normal <- p2g.normal[p2g.normal$peakType == "Distal",]
p2g.normal <- p2g.normal[p2g.normal$peakName %ni% p2g.cancer$peakName,]

df.cancer <- data.frame(num.genes = table(p2g.cancer$peakName))
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq < 3,"1-2",df.cancer$num.genes.Freq)
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq >2,"3+",df.cancer$cat)

df.normal <- data.frame(num.genes = table(p2g.normal$peakName))
df.normal$cat <- ifelse(df.normal$num.genes.Freq < 3,"1-2",df.normal$num.genes.Freq)
df.normal$cat <- ifelse(df.normal$num.genes.Freq >2,"3+",df.normal$cat)

head(df.normal)
head(df.cancer)

df.normal$type <- "normal"
df.cancer$type <- "cancer"

df<- rbind(df.normal,df.cancer)


p1 <- ggplot(df, aes(x=type, y=num.genes.Freq, fill = type)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(1.8,2.20))+NoLegend()

df.cancer %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.cancer)*100)) -> cancer
df.normal %>% 
  dplyr::count(cat) %>% 
  dplyr::mutate(perc = (n / nrow(df.normal)*100)) -> normal

cancer$type <- "cancer"
normal$type <- "normal"

comb <- rbind(cancer,normal)
print(comb)
p2 <- ggplot(comb,aes(x=cat,y=perc,fill=type))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=round(perc,3)), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+NoLegend()


CombinePlots(list(p2,p1),ncol=2)+ggsave("Barcharts_P2G_plots-HGSOC.pdf",width=6,height=3)



# Cancer specific peaks link to more genes on average with statistical significance:
test <- wilcox.test(df.cancer$num.genes.Freq,df.normal$num.genes.Freq,correct = T)
print(test)
print(test$p.value)
