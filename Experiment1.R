#######################################################################
#Experiment 1:
#Day 11 corresponds to 1DPM
#Day 16 corresponds to 4DPM

#Array Quality Metrics
###Loading packages and setting up
source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#BiocManager::install("package_name") #use to install whatever packages you might be missing

library(plyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(DESeq)
library(cowplot)
library(readr)
library(RColorBrewer)
library(gplots)
library(knitr)
library(plotly)
library(vegan)
library(kableExtra)
library(reshape2)
library(prettydoc)
library(VennDiagram)
library(MASS)
library(ggrepel)
library(stringr)
library(dplyr) #packageVersion("dplyr"); citation("dplyr")
library(pheatmap) #packageVersion("pheatmap"); citation("pheatmap")
library(adegenet) #packageVersion("adegenet"); citation("adegenet")
library(WGCNA) #packageVersion("WGCNA"); citation("WGCNA")


#conduct array quality metrics to get rid of outliers
#read in counts 
countData <- read.table("B33_counts_newtranscriptome.txt")
head(countData)
length(countData[,1])
#27008
#row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

treat=c( "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6",  "pH7.6", "pH7.6", "pH8", "pH8",  "pH8", "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6", "pH7.6", "pH8",  "pH8", "pH8","pH8", "pH7.5", "pH7.5","pH7.6", "pH7.6", "pH8", "pH8", 
         "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6","pH7.6", "pH8", "pH8", "pH8")
day=as.factor(c( "11", "11", "11", "11", "11",  "11", "11",  "11", "11", "11", "16", "16", "16", "16", "16", "16", "16", "16",  "16", "16",
                 "4", "4", "4","4", "4", "4", "8", "8", "8","8", "8", "8", "8", "8", "8"))
conditions=data.frame(treat, day)
nrow(conditions) 
ncol(countData)
real=newCountDataSet(countData,conditions) 
real=estimateSizeFactors(real)
plot(sort(sizeFactors(real))) 

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

arrayQualityMetrics(vsdBlind,intgroup=c("treat"), force=TRUE)

#6 outliers 
#Failed 2 or 3 tests: 11_7.5_3, 11_7.6_3, 16_8_1, 8_7.6_3, 8_7.6_4
#Failed 1 test: 4_7.5_2 (not super offensive, keep it in)

########################################################################
#Unloads DESeq before loading DESeq2 since the packages conflict
detach("package:DESeq", unload=TRUE)
#Load additional packages and filter data
library(DESeq2) #packageVersion("DESeq2"); citation("DESeq2")

#read in counts (old data pre-filtered)
countData <- read.table("B33_counts_newtranscriptome.txt")
head(countData)
length(countData[,1])
#27008

names(countData)=sub(".fastq.trim.sam.counts","",names(countData))
names(countData)

#removing outlier samples
countData$b33_11_7.5_3 <- NULL
countData$b33_11_7.6_3 <- NULL
countData$b33_16_8_1 <- NULL
countData$b33_8_7.6_3 <- NULL
countData$b33_8_7.6_4 <- NULL

totalCounts=colSums(countData)
totalCounts
barplot(totalCounts)
# b33_11_7.5_1 b33_11_7.5_2 b33_11_7.6_1 b33_11_7.6_2 b33_11_7.6_4   b33_11_8_1 
# 972423       750239      1156530      1008800       905133      1404888 
# b33_11_8_2   b33_11_8_4 b33_16_7.5_2 b33_16_7.5_3 b33_16_7.5_4 b33_16_7.6_1 
# 1154117       887547       733385      1459340      1255646      1299882 
# b33_16_7.6_3 b33_16_7.6_4   b33_16_8_2   b33_16_8_3   b33_16_8_4  b33_4_7.5_1 
# 1146724      1557855      1070388       987887      1381910      1142969 
# b33_4_7.5_2  b33_4_7.6_2  b33_4_7.6_3    b33_4_8_2    b33_4_8_3  b33_8_7.5_2 
# 665524      1327316       932421       954244      1468909      1003391 
# b33_8_7.5_3  b33_8_7.5_4  b33_8_7.6_2   b33_8_8_12    b33_8_8_3    b33_8_8_4 
# 1368268       918446       330871       902746      1055424      1063168 

min(totalCounts) #330,871
max(totalCounts)  # 1,557,855

treat=c( "pH7.5", "pH7.5", "pH7.6", "pH7.6",  "pH7.6", "pH8", "pH8",  "pH8", "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6", "pH7.6", "pH8",  "pH8", "pH8", 
         "pH7.5", "pH7.5","pH7.6", "pH7.6", "pH8", "pH8", 
         "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH8", "pH8", "pH8")
day=as.factor(c( "11", "11", "11", "11", "11",  "11", "11",  "11", "16", "16", "16", "16", "16", "16", "16", "16",  "16",
                 "4", "4", "4","4", "4", "4",  
                 "8", "8", "8","8", "8", "8", "8"))

#Dividing data up into time points
colData =data.frame(treat, day)
colData4 = colData[colData$day == 4,]
colData8 = colData[colData$day == 8,]
colDatalarvae = rbind(colData4, colData8)
colData11 = colData[colData$day == 11,]
colData16 = colData[colData$day == 16,]
colDatajuv = rbind(colData11, colData16)

#Creating a Deseq matrix
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat+day+treat*day ) 
head(dds)
dds4 = dds[,colData$day == 4]
head(dds4)
dds8 = dds[,colData$day == 8]
head(dds8)
ddslarvae = dds[,colData$day == 4 | colData$day == 8] 
head(ddslarvae)
dds11 = dds[,colData$day == 11]
head(dds11)
dds16 = dds[,colData$day == 16]
head(dds16)
ddsjuvs = dds[,colData$day == 11 | colData$day == 16] 
head(ddsjuvs)

design(dds4) = ~treat 
design(dds8) = ~treat
design(ddslarvae) = ~treat
design(dds11) = ~treat
design(dds16) = ~treat
design(ddsjuvs) = ~treat

#Running matrices through DESeq
dds <-DESeq(dds)
dds4<-DESeq(dds4)
dds8<-DESeq(dds8)
dds11<-DESeq(dds11)
dds16<-DESeq(dds16)
ddslarvae<-DESeq(ddslarvae)
ddsjuvs<-DESeq(ddsjuvs)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot")
plotDispEsts(dds4, main="Dispersion plot (Day 4)")
plotDispEsts(dds8, main="Dispersion plot (Day 8)")
plotDispEsts(dds11, main="Dispersion plot (Day 11)")
plotDispEsts(dds16, main="Dispersion plot (Day 16)")

#Rlog results
dds4logged = DESeq2::rlog(dds4, blind = TRUE)
dds8logged = DESeq2::rlog(dds8, blind = TRUE)
dds11logged = DESeq2::rlog(dds11, blind = TRUE)
dds16logged = DESeq2::rlog(dds16, blind = TRUE)
ddslarvaelogged = DESeq2::rlog(ddslarvae, blind = TRUE)
ddsjuvslogged = DESeq2::rlog(ddsjuvs, blind = TRUE)

#Saving Rlog results
write.csv(assay(dds4logged), "dds4_rlogged.csv")
write.csv(assay(dds8logged), "dds8_rlogged.csv")
write.csv(assay(dds11logged), "dds11_rlogged.csv")
write.csv(assay(dds16logged), "dds16_rlogged.csv")
write.csv(assay(ddslarvaelogged), "ddslarvae_rlogged.csv")
write.csv(assay(ddsjuvslogged), "ddsjuvs_rlogged.csv")

#Separating samples into individual columns
dds4rlogged = read.table("dds4_rlogged.csv")
dds4rlogged = as.data.frame(dds4rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7"), ",")
  dds4rlogged$`1` = NULL
  rownames(dds4rlogged) = dds4rlogged$V1
  colnames(dds4rlogged) = dds4rlogged[1,]
  colnames(dds4rlogged)[1] = "Gene"
 

dds8rlogged = read.table("dds8_rlogged.csv")
dds8rlogged = as.data.frame(dds8rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7", "8"), ",")
dds8rlogged$`1` = NULL
rownames(dds8rlogged) = dds8rlogged$V1
colnames(dds8rlogged) = dds8rlogged[1,]
colnames(dds8rlogged)[1] = "Gene"
  
  
dds11rlogged = read.table("dds11_rlogged.csv")
dds11rlogged = as.data.frame(dds11rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7", "8","9"), ",")
dds11rlogged$`1` = NULL
rownames(dds11rlogged) = dds11rlogged$V
colnames(dds11rlogged) = dds11rlogged[1,]
colnames(dds11rlogged)[1] = "Gene"

dds16rlogged = read.table("dds16_rlogged.csv")
dds16rlogged = as.data.frame(dds16rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7","8","9", "10"), ",")
dds16rlogged$`1` = NULL
rownames(dds16rlogged) = dds16rlogged$V
colnames(dds16rlogged) = dds16rlogged[1,]
colnames(dds16rlogged)[1] = "Gene"

ddslarvaerlogged = read.table("ddslarvae_rlogged.csv")
ddsjuvsrlogged = read.table("ddsjuvs_rlogged.csv")
#the same code can be applied to the larvae and juvs subsets to
#separate the rlog results into columns. 

#DESEq results w/ and w/o call to contrast by pH
#also merging DESEq and rlog results into one table
res4<- results(dds4)
res4 = as.data.frame(res4)
res4 = merge(res4, dds4rlogged, by=0)
rownames(res4) <- res4$Gene
res4_7.5<- results(dds4, contrast=c("treat","pH7.5","pH8"))
res4_7.5 = as.data.frame(res4_7.5)
res4_7.5 = merge(res4_7.5, dds4rlogged, by=0)
rownames(res4_7.5) <- res4_7.5$Gene
res4_7.6<- results(dds4, contrast=c("treat","pH7.6","pH8"))
res4_7.6 = as.data.frame(res4_7.6)
res4_7.6 = merge(res4_7.6, dds4rlogged, by=0)
rownames(res4_7.6) <- res4_7.6$Gene

res8<- results(dds8)
res8 = as.data.frame(res8)
res8 = merge(res8, dds8rlogged, by=0)
rownames(res8) <- res8$Gene
res8_7.5<- results(dds8, contrast=c("treat","pH7.5","pH8"))
res8_7.5 = as.data.frame(res8_7.5)
res8_7.5 = merge(res8_7.5, dds8rlogged, by=0)
rownames(res8_7.5) <- res8_7.5$Gene
res8_7.6<- results(dds8, contrast=c("treat","pH7.6","pH8"))
res8_7.6 = as.data.frame(res8_7.6)
res8_7.6 = merge(res8_7.6, dds8rlogged, by=0)
rownames(res8_7.6) <- res8_7.6$Gene

res11<- results(dds11)
res11 = as.data.frame(res11)
res11 = merge(res11, dds11rlogged, by=0)
rownames(res11) <- res11$Gene
res11_7.5<- results(dds11, contrast=c("treat","pH7.5","pH8"))
res11_7.5 = as.data.frame(res11_7.5)
res11_7.5 = merge(res11_7.5, dds11rlogged, by=0)
rownames(res11_7.5) <- res11_7.5$Gene
res11_7.6<- results(dds11, contrast=c("treat","pH7.6","pH8"))
res11_7.6 = as.data.frame(res11_7.6)
res11_7.6 = merge(res11_7.6, dds11rlogged, by=0)
rownames(res11_7.6) <- res11_7.6$Gene

res16<- results(dds16)
res16 = as.data.frame(res16)
res16 = merge(res16, dds16rlogged, by=0)
rownames(res16) <- res16$Gene
res16_7.5<- results(dds16, contrast=c("treat","pH7.5","pH8"))
res16_7.5 = as.data.frame(res16_7.5)
res16_7.5 = merge(res16_7.5, dds16rlogged, by=0)
rownames(res16_7.5) <- res16_7.5$Gene
res16_7.6<- results(dds16, contrast=c("treat","pH7.6","pH8"))
res16_7.6 = as.data.frame(res16_7.6)
res16_7.6 = merge(res16_7.6, dds16rlogged, by=0)
rownames(res16_7.6) <- res16_7.6$Gene

reslarvae = as.data.frame(results(ddslarvae))
reslarvae = merge(reslarvae, ddslarvaerlogged, by=0)
rownames(reslarvae) <- reslarvae$Gene
write.csv(reslarvae, "Larvae_Results.csv")

resjuvs = as.data.frame(results(ddsjuvs))
resjuvs = merge(resjuvs, ddsjuvsrlogged, by=0)
rownames(resjuvs) <- resjuvs$Gene
write.csv(resjuvs, "Juvs_Results.csv")

write.csv(res16, "Day16_Results.csv", row.names = TRUE)
write.csv(res4_7.6, "Day4_7.6_Results.csv", row.names = TRUE)
write.csv(res8_7.5, "Day8_7.5_Results.csv", row.names = TRUE)
write.csv(res8_7.6, "Day8_7.6_Results.csv", row.names = TRUE)
write.csv(res11_7.5, "Day11_7.5_Results.csv", row.names = TRUE)
write.csv(res11_7.6, "Day11_7.6_Results.csv", row.names = TRUE)
write.csv(res16_7.5, "Day16_7.5_Results.csv", row.names = TRUE)
write.csv(res16_7.6, "Day16_7.6_Results.csv", row.names = TRUE)

#############################################################################
#############################################################################
#Pairwise comparisons/# of DEGs: pH 7.6 treatments
respH76.4 <- res4_7.6
respH76.8 <- res8_7.6
respH76.11 <- res11_7.6
respH76.16 <- res16_7.6

table(respH76.4$padj<.1)
table(respH76.4$padj<.05)
table(respH76.4$padj<.01)
# 0.1=0/26651
# 0.05=0/26651
# 0.01=0/26651
summary(respH76.4)
#LFC>0: 0, 0%
#LFC<0: 0, 0%
#Outliers: 0, 0%
#Low counts: 0, 0%

table(respH76.8$padj<.1)
table(respH76.8$padj<.05)
table(respH76.8$padj<.01)
# 0.1=0/26914
# 0.05=0/26914
# 0.01=0/26914
summary(respH76.8)
#LFC>0: 0, 0%
#LFC<0: 0, 0%
#Outliers: 6, 0.022%
#Low counts: 0, 0%

table(respH76.11$padj<.1)
table(respH76.11$padj<.05)
table(respH76.11$padj<.01)
# 0.1=23/13320
# 0.05=17/13320
# 0.01=7/13320
summary(respH76.11)
#LFC>0: 2, 0.0074%
#LFC<0: 21, 0.0078%
#Outliers: 49, 0.18%
#Low counts: 13582, 50%

table(respH76.16$padj<.1)
table(respH76.16$padj<.05)
table(respH76.16$padj<.01)
# 0.1=1/26957
# 0.05=1/26957
# 0.01=0/26957
summary(respH76.16)
#LFC>0: 1, 0.0037%
#LFC<0: 0, 0%
#Outliers: 26, 0.096%
#Low counts: 0, 0%

nrow(respH76.4[respH76.4$padj<0.05 & !is.na(respH76.4$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0
nrow(respH76.8[respH76.8$padj<0.05 & !is.na(respH76.8$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0
nrow(respH76.11[respH76.11$padj<0.05 & !is.na(respH76.11$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #17
nrow(respH76.16[respH76.16$padj<0.05 & !is.na(respH76.16$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #1

plotMA(respH76.4, main="pH8 vs pH7.6 (Day 4)", ylim=c(-1,1))
plotMA(respH76.8, main="pH8 vs pH7.6 (Day 8)", ylim=c(-1,1))
plotMA(respH76.11, main="pH8 vs pH7.6 (1 DPM)", ylim=c(-1,1))
plotMA(respH76.16, main="pH8 vs pH7.6 (4 DPM)", ylim=c(-1,1))

nrow(respH76.4[respH76.4$padj<0.1 & respH76.4$log2FoldChange > 0 & !is.na(respH76.4$padj),])
nrow(respH76.4[respH76.4$padj<0.1 & respH76.4$log2FoldChange < -0 & !is.na(respH76.4$padj),])
#UP in 7.6 (day 4): 0
#DOWN in 7.6 (day 4): 0
nrow(respH76.8[respH76.8$padj<0.1 & respH76.8$log2FoldChange > 0 & !is.na(respH76.8$padj),])
nrow(respH76.8[respH76.8$padj<0.1 & respH76.8$log2FoldChange < -0 & !is.na(respH76.8$padj),])
#UP in 7.6 (day 8): 0
#DOWN in 7.6 (day 8): 0
nrow(respH76.11[respH76.11$padj<0.1 & respH76.11$log2FoldChange > 0 & !is.na(respH76.11$padj),])
nrow(respH76.11[respH76.11$padj<0.1 & respH76.11$log2FoldChange < -0 & !is.na(respH76.11$padj),])
#UP in 7.6 (day 11): 2
#DOWN in 7.6 (day 11): 21
nrow(respH76.16[respH76.16$padj<0.1 & respH76.16$log2FoldChange > 0 & !is.na(respH76.16$padj),])
nrow(respH76.16[respH76.16$padj<0.1 & respH76.16$log2FoldChange < -0 & !is.na(respH76.16$padj),])
#UP in 7.6 (day 16): 1
#DOWN in 7.6 (day 16): 0

write.table(respH76.4, file="7.6.4_2017.txt", quote=F, sep="\t")
write.table(respH76.8, file="7.6.8_2017.txt", quote=F, sep="\t")
write.table(respH76.11, file="7.6.11_2017.txt", quote=F, sep="\t")
write.table(respH76.16, file="7.6.16_2017.txt", quote=F, sep="\t")

MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}

#orange dots = upregulated, blue dots = downregulated
MAttPlot(respH76.4)
MAttPlot(respH76.8)
MAttPlot(respH76.11)
MAttPlot(respH76.16)

cd4 <- read.table("7.6.4_2017.txt")
head(cd4)
cd8 <- read.table("7.6.8_2017.txt")
head(cd8)
cd11 <- read.table("7.6.11_2017.txt")
head(cd11)
cd16 <- read.table("7.6.16_2017.txt")
head(cd16)

##make the GO table for MWU
cd4$isogroup=row.names(cd4)
cd8$isogroup=row.names(cd8)
cd11$isogroup=row.names(cd11)
cd16$isogroup=row.names(cd16)

go_input_7.6.4 = cd4 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.6.4)
colnames(go_input_7.6.4) <- c("gene", "pval")
head(go_input_7.6.4)
write.csv(go_input_7.6.4, file="7.6.4_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.6.8 = cd8 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.6.8)
colnames(go_input_7.6.8) <- c("gene", "pval")
head(go_input_7.6.8)
write.csv(go_input_7.6.8, file="7.6.8_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.6.11 = cd11 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.6.11)
colnames(go_input_7.6.11) <- c("gene", "pval")
head(go_input_7.6.11)
write.csv(go_input_7.6.11, file="7.6.11_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.6.16 = cd16 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.6.16)
colnames(go_input_7.6.16) <- c("gene", "pval")
head(go_input_7.6.16)
write.csv(go_input_7.6.16, file="7.6.16_rec_GO.csv", quote=F, row.names=FALSE)
##################################################################################################
##################################################################################################
#Pairwise comparisons/# of DEGs: pH 7.5 treatments
respH75.4 <- res4_7.5
respH75.8 <- res8_7.5
respH75.11 <- res11_7.5
respH75.16 <- res16_7.5

table(respH75.4$padj<.1)
table(respH75.4$padj<.05)
table(respH75.4$padj<.01)
# 0.1=0/26651
# 0.05=0/26651
# 0.01=0/26651
summary(respH75.4)
#LFC>0: 0, 0%
#LFC<0: 0, 0%
#Outliers: 0, 0%
#Low counts: 0, 0%

table(respH75.8$padj<.1)
table(respH75.8$padj<.05)
table(respH75.8$padj<.01)
# 0.1=1/26914
# 0.05=1/26914
# 0.01=1/26914
summary(respH75.8)
#LFC>0: 0, 0%
#LFC<0: 1, 0.0037%
#Outliers: 6, 0.022%
#Low counts: 0, 0%

table(respH75.11$padj<.1)
table(respH75.11$padj<.05)
table(respH75.11$padj<.01)
# 0.1=92/11756
# 0.05=65/11756
# 0.01=30/11756
summary(respH75.11)
#LFC>0: 65, 0.24%
#LFC<0: 27, 0.1%
#Outliers: 49, 0.18%
#Low counts: 15146, 56%

table(respH75.16$padj<.1)
table(respH75.16$padj<.05)
table(respH75.16$padj<.01)
# 0.1=210/18064
# 0.05=121/18064
# 0.01=51/18064
summary(respH75.16)
#LFC>0: 140, 0.52%
#LFC<0: 70, 0.26%
#Outliers: 26, 0.096%
#Low counts: 8893, 33%

nrow(respH75.4[respH75.4$padj<0.05 & !is.na(respH75.4$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0
nrow(respH75.8[respH75.8$padj<0.05 & !is.na(respH75.8$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #1
nrow(respH75.11[respH75.11$padj<0.05 & !is.na(respH75.11$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #65
nrow(respH75.16[respH75.16$padj<0.05 & !is.na(respH75.16$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #121

plotMA(respH75.4, main="pH8 vs pH7.5 (Day 4)", ylim=c(-1,1))
plotMA(respH75.8, main="pH8 vs pH7.5 (Day 8)", ylim=c(-1,1))
plotMA(respH75.11, main="pH8 vs pH7.5 (1 DPM)", ylim=c(-1,1))
plotMA(respH75.16, main="pH8 vs pH7.5 (4 DPM)", ylim=c(-1,1))

nrow(respH75.4[respH75.4$padj<0.1 & respH75.4$log2FoldChange > 0 & !is.na(respH75.4$padj),])
nrow(respH75.4[respH75.4$padj<0.1 & respH75.4$log2FoldChange < -0 & !is.na(respH75.4$padj),])
#UP in 7.5 (day 4): 0
#DOWN in 7.5 (day 4): 0
nrow(respH75.8[respH75.8$padj<0.1 & respH75.8$log2FoldChange > 0 & !is.na(respH75.8$padj),])
nrow(respH75.8[respH75.8$padj<0.1 & respH75.8$log2FoldChange < -0 & !is.na(respH75.8$padj),])
#UP in 7.5 (day 8): 0
#DOWN in 7.5 (day 8): 1
nrow(respH75.11[respH75.11$padj<0.1 & respH75.11$log2FoldChange > 0 & !is.na(respH75.11$padj),])
nrow(respH75.11[respH75.11$padj<0.1 & respH75.11$log2FoldChange < -0 & !is.na(respH75.11$padj),])
#UP in 7.5 (day 11): 65
#DOWN in 7.5 (day 11): 27
nrow(respH75.16[respH75.16$padj<0.1 & respH75.16$log2FoldChange > 0 & !is.na(respH75.16$padj),])
nrow(respH75.16[respH75.16$padj<0.1 & respH75.16$log2FoldChange < -0 & !is.na(respH75.16$padj),])
#UP in 7.5 (day 16): 140
#DOWN in 7.5 (day 16): 70

write.table(respH75.4, file="7.5.4_2017.txt", quote=F, sep="\t")
write.table(respH75.8, file="7.5.8_2017.txt", quote=F, sep="\t")
write.table(respH75.11, file="7.5.11_2017.txt", quote=F, sep="\t")
write.table(respH75.16, file="7.5.16_2017.txt", quote=F, sep="\t")

cp4 <- read.table("7.5.4_2017.txt")
head(cp4)
cp8 <- read.table("7.5.8_2017.txt")
head(cp8)
cp11 <- read.table("7.5.11_2017.txt")
head(cp11)
cp16 <- read.table("7.5.16_2017.txt")
head(cp16)

MAttPlot(cp4)
MAttPlot(cp8)
MAttPlot(cp11)
MAttPlot(cp16)

##make the GO table for MWU
cp4$isogroup=row.names(cp4)
cp8$isogroup=row.names(cp8)
cp11$isogroup=row.names(cp11)
cp16$isogroup=row.names(cp16)

go_input_7.5.4 = cp4 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.4)
colnames(go_input_7.5.4) <- c("gene", "pval")
head(go_input_7.5.4)
write.csv(go_input_7.5.4, file="7.5.4_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.5.8 = cp8 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.8)
colnames(go_input_7.5.8) <- c("gene", "pval")
head(go_input_7.5.8)
write.csv(go_input_7.5.8, file="7.5.8_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.5.11 = cp11 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.11)
colnames(go_input_7.5.11) <- c("gene", "pval")
head(go_input_7.5.11)
write.csv(go_input_7.5.11, file="7.5.11_rec_GO.csv", quote=F, row.names=FALSE)

go_input_7.5.16 = cp16 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.16)
colnames(go_input_7.5.16) <- c("gene", "pval")
head(go_input_7.5.16)
write.csv(go_input_7.5.16, file="7.5.16_rec_GO.csv", quote=F, row.names=FALSE)

################################################################################
#PCAs
#This is the code used to make the PCA for juveniles overall. It can be used as a template
#for making all the other PCAs
PCA_colors <- ifelse(colDatajuv$treat == 'pH7.5', 'darkgoldenrod2', 
                     ifelse(colDatajuv$treat == 'pH8','dodgerblue3', 'darkolivegreen3'))
Pt_colors = c("darkgoldenrod2","darkolivegreen3", "dodgerblue3")
juv_colors <- ifelse(colDatajuv$day == '11', 'turquoise1',
                     ifelse(colDatajuv$day == '16', 'slateblue1','black'))
pt_sym <- ifelse(colDatajuv$treat == 'pH7.5', 15, 
                 ifelse(colDatajuv$treat == 'pH8',17, 16))
juvs_pcadata = DESeq2::plotPCA(ddsjuvslogged, intgroup ="treat", returnData = TRUE)
juvs_percentVar = round(100* attr(juvs_pcadata, "percentVar"))
juvs_pca =prcomp(t(assay(ddsjuvslogged)), center = TRUE, scale. = FALSE)

PCA_juvs =as.data.frame(juvs_pca$x)%>%
  dplyr::select(PC1, PC2)


plot(PCA_juvs,choices=c(1,2),display="sites",type="n")
points(PCA_juvs,col=juv_colors, pch=pt_sym)
#ordispider(ccaa,groups=treat,col=c("green", "blue"))
#ordispider(ccaa,groups=conds$hour,col=c("green", "blue", "red", "yellow"))
ordihull(PCA_juvs,groups=colDatajuv$treat,draw="polygon",col="grey95",label=T)
ordihull(PCA_juvs,groups=colDatajuv$treat,draw="polygon",border=Pt_colors,label=F)
legend(x="topright", c("1DPM", "4DPM"), fill= c("turquoise1","slateblue1"))

adonis(juvs_pca$x ~ colDatajuv$treat + colDatajuv$day, method='eu')

###################################################################################
###########GO Enrichment Analysis
#GO enrichment analysis

#Day 4
#pH 7.5 v. 8

#######MF
# Edit these to match your data file names: 
input="7.5.4_rec_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="crepidula_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#pH 7.6 v 8

#######MF
input="7.6.4_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#Day 8
#pH 7.5 v 8

#######MF
input="7.5.8_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#pH 7.6 v 8

#######MF
input="7.6.8_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#Day 11
#pH 7.5 v 8

#######MF
input="7.5.11_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#pH 7.6 v 8

#######MF
input="7.6.11_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#Day 16
#pH 7.5 v 8

#######MF
input="7.5.16_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#pH 7.6 v 8

#######MF
input="7.6.16_rec_GO.csv" 
goDivision="MF" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#######BP
goDivision="BP" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######CC
goDivision="CC" # either MF, or BP, or CC

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25
)

# Plotting results
results=gomwuPlot(input,goAnnotations,goDivision,
                  absValue=1,
                  level1=0.1,
                  level2=0.05,
                  level3=0.01,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

##############################################################################
#####################################################################
#DEG heatmaps
library(tidyverse)
iso2gene = read.delim("crepidula_iso2gene.tab", header=FALSE)
rownames(iso2gene) <- iso2gene$V1


#Day 11 - 7.5
D11_7.5_Results = read.csv("Day11_7.5_Results.csv")
rownames(D11_7.5_Results) = D11_7.5_Results$X
D11_7.5_Results$X = NULL
rld_data_D11_7.5=as.data.frame(D11_7.5_Results[,c((5:6),(8:9), (13:15))])

p.val=0.10 # FDR cutoff
conds_D11_7.5=rld_data_D11_7.5[rld_data_D11_7.5$padj<=p.val & !is.na(rld_data_D11_7.5$padj),]
length(conds_D11_7.5[,1])
#92
exp_D11_7.5=conds_D11_7.5[,3:7] #change to your columns
head(exp_D11_7.5)
means_D11_7.5=apply(exp_D11_7.5,1,mean) # means of rows
explc_D11_7.5=exp_D11_7.5-means_D11_7.5 # subtracting them


D11_7.5_df_all_iso <- explc_D11_7.5 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))

D11_7.5_unanno=D11_7.5_df_all_iso[,2:7]

D11_7.5_df_only_anno <- D11_7.5_df_all_iso %>%
  filter(!is.na(V2))
rownames(D11_7.5_df_only_anno) <- make.unique(D11_7.5_df_only_anno$V2)
D11_7.5_df_only_anno$V1 = NULL

D11_7.5_df_unanno <- D11_7.5_df_all_iso %>%
  filter(is.na(V2))
rownames(D11_7.5_df_unanno) <- D11_7.5_df_unanno$V1
D11_7.5_df_unanno$V1 = NULL
D11_7.5_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D11_7.5_df_only_anno)

D11_7.5_my_sample_col <- data.frame(Source = c("pH 7.5", "pH 7.5", "pH 8", "pH 8", "pH 8"))
D11_7.5_df_only_anno$V2 = NULL
colnames(D11_7.5_df_only_anno) = row.names(D11_7.5_my_sample_col)


my_colour_7.5 = list(
  Source = c('pH 7.5' = 'darkgoldenrod2', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D11_7.5_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D11_7.5_my_sample_col, annotation_colors =my_colour_7.5, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D11_7.5_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D11_7.5_my_sample_col, annotation_colors =my_colour_7.5, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

#Day 11 - 7.6
D11_7.6_Results = read.csv("Day11_7.6_Results.csv")
rownames(D11_7.6_Results) = D11_7.6_Results$X
D11_7.6_Results$X = NULL
rld_data_D11_7.6=as.data.frame(D11_7.6_Results[,c((5:6),(10:15))])
p.val=0.10 # FDR cutoff
conds_D11_7.6=rld_data_D11_7.6[rld_data_D11_7.6$padj<=p.val & !is.na(rld_data_D11_7.6$padj),]
length(conds_D11_7.6[,1])
#23
exp_D11_7.6=conds_D11_7.6[,3:8] #change to your columns
head(exp_D11_7.6)
means_D11_7.6=apply(exp_D11_7.6,1,mean) # means of rows
explc_D11_7.6=exp_D11_7.6-means_D11_7.6 # subtracting them


D11_7.6_df_all_iso <- explc_D11_7.6 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))


D11_7.6_df_only_anno <- D11_7.6_df_all_iso %>%
  filter(!is.na(V2))
rownames(D11_7.6_df_only_anno) <- make.unique(D11_7.6_df_only_anno$V2)
D11_7.6_df_only_anno$V1 = NULL

D11_7.6_df_unanno <- D11_7.6_df_all_iso %>%
  filter(is.na(V2))
rownames(D11_7.6_df_unanno) <- D11_7.6_df_unanno$V1
D11_7.6_df_unanno$V1 = NULL
D11_7.6_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D11_7.6_df_only_anno)

D11_7.6_my_sample_col <- data.frame(Source = c("pH 7.6", "pH 7.6","pH 7.6", "pH 8", "pH 8", "pH 8"))
D11_7.6_df_only_anno$V2 = NULL
row.names(D11_7.6_my_sample_col) = colnames(D11_7.6_df_only_anno) 


my_colour_7.6 = list(
  Source = c('pH 7.6' = 'darkolivegreen3', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D11_7.6_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D11_7.6_my_sample_col, annotation_colors =my_colour_7.6, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D11_7.6_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D11_7.6_my_sample_col, annotation_colors =my_colour_7.6, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

#Day 11 - overall
D11_Results = read.csv("Day11_Results.csv")
rownames(D11_Results) = D11_Results$X
D11_Results$X = NULL
rld_data_D11=as.data.frame(D11_Results[,c((6:7),(10:17))])
p.val=0.10 # FDR cutoff
conds_D11=rld_data_D11[rld_data_D11$padj<=p.val & !is.na(rld_data_D11$padj),]
length(conds_D11[,1])
#92
exp_D11=conds_D11[,3:10] #change to your columns
head(exp_D11)
means_D11=apply(exp_D11,1,mean) # means of rows
explc_D11=exp_D11-means_D11 # subtracting them


D11_df_all_iso <- explc_D11 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))


D11_df_only_anno <- D11_df_all_iso %>%
  filter(!is.na(V2))
rownames(D11_df_only_anno) <- make.unique(D11_df_only_anno$V2)
D11_df_only_anno$V1 = NULL

D11_df_unanno <- D11_df_all_iso %>%
  filter(is.na(V2))
rownames(D11_df_unanno) <- D11_df_unanno$V1
D11_df_unanno$V1 = NULL
D11_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D11_df_only_anno)

D11_my_sample_col <- data.frame(Source = c("pH 7.5","pH 7.5","pH 7.6", "pH 7.6","pH 7.6", "pH 8", "pH 8", "pH 8"))
D11_df_only_anno$V2 = NULL
row.names(D11_my_sample_col) = colnames(D11_df_only_anno) 


my_colour_D11 = list(
  Source = c('pH 7.5' = 'darkgoldenrod2','pH 7.6' = 'darkolivegreen3', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D11_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D11_my_sample_col, annotation_colors =my_colour_D11, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D11_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D11_my_sample_col, annotation_colors =my_colour_D11, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

#Day 16 - 7.5
D16_7.5_Results = read.csv("Day16_7.5_Results.csv")
rownames(D16_7.5_Results) = D16_7.5_Results$X
D16_7.5_Results$X = NULL
rld_data_D16_7.5=as.data.frame(D16_7.5_Results[,c((5:6),(8:10), (14:16))])

p.val=0.10 # FDR cutoff
conds_D16_7.5=rld_data_D16_7.5[rld_data_D16_7.5$padj<=p.val & !is.na(rld_data_D16_7.5$padj),]
length(conds_D16_7.5[,1])
#210
exp_D16_7.5=conds_D16_7.5[,3:8] #change to your columns
head(exp_D16_7.5)
means_D16_7.5=apply(exp_D16_7.5,1,mean) # means of rows
explc_D16_7.5=exp_D16_7.5-means_D16_7.5 # subtracting them


D16_7.5_df_all_iso <- explc_D16_7.5 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))

D16_7.5_unanno=D16_7.5_df_all_iso[,2:8]

D16_7.5_df_only_anno <- D16_7.5_df_all_iso %>%
  filter(!is.na(V2))
rownames(D16_7.5_df_only_anno) <- make.unique(D16_7.5_df_only_anno$V2)
D16_7.5_df_only_anno$V1 = NULL

D16_7.5_df_unanno <- D16_7.5_df_all_iso %>%
  filter(is.na(V2))
rownames(D16_7.5_df_unanno) <- D16_7.5_df_unanno$V1
D16_7.5_df_unanno$V1 = NULL
D16_7.5_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D16_7.5_df_only_anno)

D16_7.5_my_sample_col <- data.frame(Source = c("pH 7.5", "pH 7.5","pH 7.5", "pH 8", "pH 8", "pH 8"))
D16_7.5_df_only_anno$V2 = NULL
colnames(D16_7.5_df_only_anno)=row.names(D16_7.5_my_sample_col)
row.names(D16_7.5_my_sample_col)=colnames(D16_7.5_df_only_anno)

my_colour_7.5 = list(
  Source = c('pH 7.5' = 'darkgoldenrod2', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D16_7.5_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D16_7.5_my_sample_col, annotation_colors =my_colour_7.5, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D16_7.5_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D16_7.5_my_sample_col, annotation_colors =my_colour_7.5, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

#Day 16 - 7.6
D16_7.6_Results = read.csv("Day16_7.6_Results.csv")
rownames(D16_7.6_Results) = D16_7.6_Results$X
D16_7.6_Results$X = NULL
rld_data_D16_7.6=as.data.frame(D16_7.6_Results[,c((5:6),(11:16))])
p.val=0.10 # FDR cutoff
conds_D16_7.6=rld_data_D16_7.6[rld_data_D16_7.6$padj<=p.val & !is.na(rld_data_D16_7.6$padj),]
length(conds_D16_7.6[,1])
#1
exp_D16_7.6=conds_D16_7.6[,3:8] #change to your columns
head(exp_D16_7.6)
means_D16_7.6=apply(exp_D16_7.6,1,mean) # means of rows
explc_D16_7.6=exp_D16_7.6-means_D16_7.6 # subtracting them


D16_7.6_df_all_iso <- explc_D16_7.6 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))


D16_7.6_df_only_anno <- D16_7.6_df_all_iso %>%
  filter(!is.na(V2))
rownames(D16_7.6_df_only_anno) <- make.unique(D16_7.6_df_only_anno$V2)
D16_7.6_df_only_anno$V1 = NULL

D16_7.6_df_unanno <- D16_7.6_df_all_iso %>%
  filter(is.na(V2))
rownames(D16_7.6_df_unanno) <- D16_7.6_df_unanno$V1
D16_7.6_df_unanno$V1 = NULL
D16_7.6_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D16_7.6_df_only_anno)

D16_7.6_my_sample_col <- data.frame(Source = c("pH 7.6", "pH 7.6","pH 7.6", "pH 8", "pH 8", "pH 8"))
D16_7.6_df_only_anno$V2 = NULL
row.names(D16_7.6_my_sample_col) = colnames(D16_7.6_df_only_anno) 


my_colour_7.6 = list(
  Source = c('pH 7.6' = 'darkolivegreen3', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D16_7.6_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D16_7.6_my_sample_col, annotation_colors =my_colour_7.6, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D16_7.6_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D16_7.6_my_sample_col, annotation_colors =my_colour_7.6, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

#Day 16 - overall
D16_Results = read.csv("Day16_Results.csv")
rownames(D16_Results) = D16_Results$X
D16_Results$X = NULL
rld_data_D16=as.data.frame(D16_Results[,c((6:7),(10:18))])
p.val=0.10 # FDR cutoff
conds_D16=rld_data_D16[rld_data_D16$padj<=p.val & !is.na(rld_data_D16$padj),]
length(conds_D16[,1])
#210
exp_D16=conds_D16[,3:11] #change to your columns
head(exp_D16)
means_D16=apply(exp_D16,1,mean) # means of rows
explc_D16=exp_D16-means_D16 # subtracting them


D16_df_all_iso <- explc_D16 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))


D16_df_only_anno <- D16_df_all_iso %>%
  filter(!is.na(V2))
rownames(D16_df_only_anno) <- make.unique(D16_df_only_anno$V2)
D16_df_only_anno$V1 = NULL

D16_df_unanno <- D16_df_all_iso %>%
  filter(is.na(V2))
rownames(D16_df_unanno) <- D16_df_unanno$V1
D16_df_unanno$V1 = NULL
D16_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(D16_df_only_anno)

D16_my_sample_col <- data.frame(Source = c("pH 7.5","pH 7.5","pH 7.5","pH 7.6", "pH 7.6","pH 7.6", "pH 8", "pH 8", "pH 8"))
D16_df_only_anno$V2 = NULL
row.names(D16_my_sample_col) = colnames(D16_df_only_anno) 


my_colour_D16 = list(
  Source = c('pH 7.5' = 'darkgoldenrod2','pH 7.6' = 'darkolivegreen3', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(D16_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= D16_my_sample_col, annotation_colors =my_colour_D16, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(D16_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= D16_my_sample_col, annotation_colors =my_colour_D16, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()