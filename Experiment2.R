#######################################################################
###############################################################
#Experiment 2

#Array Quality Metrics
###Loading packages and setting up
source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#BiocManager::install("affycoretools", dependencies=TRUE) #use to install whatever packages you might be missing

library(DESeq); #packageVersion("DESeq"); citation("DESeq")
library(affycoretools); #packageVersion("affycoretools"); citation("affycoretools")
library(arrayQualityMetrics); #packageVersion("arrayQualityMetrics"); citation("arrayQualityMetrics")
library(genefilter); #packageVersion("genefilter"); citation("genefilter")

#conduct array quality metrics to get rid of outliers
#read in counts 
countData <- read.table("B37_counts_newtranscriptome.txt")
head(countData)
length(countData[,1])
#20788
#row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

treat=c( "pH8", "pH8",  "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8",  "pH8", "pH8")
hour=as.factor(c( "0", "0", "0", "0", "10",  "10", "10",  "10", "10", "10", "10", "10", "24", "24", "24", "24", "24", "24",  "24", "24",
                  "48", "48", "48","48", "48", "48", "48", "48", "4","4", "4", "4", "4", "4", "4", "4"))
conditions=data.frame(treat, hour)
nrow(conditions) 
ncol(countData)
real=newCountDataSet(countData,conditions) 
real=estimateSizeFactors(real)
plot(sort(sizeFactors(real))) 

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

arrayQualityMetrics(vsdBlind,intgroup=c("treat"), force=TRUE)

#Unloads DESeq before loading DESeq2 since the packages conflict
detach("package:DESeq", unload=TRUE)
#Load additional packages and filter data
library(DESeq2) #packageVersion("DESeq2"); citation("DESeq2")
library(ggplot2) #packageVersion("ggplot2"); citation("ggplot2")
library(dplyr) #packageVersion("dplyr"); citation("dplyr")
library(RColorBrewer) #packageVersion("RColorBrewer"); citation("RColorBrewer")
library(gplots) #packageVersion("gplots"); citation("gplots")
library(pheatmap) #packageVersion("pheatmap"); citation("pheatmap")
library(vegan) #packageVersion("vegan"); citation("vegan")
library(ggrepel) #packageVersion("ggrepel"); citation("ggrepel")
library(tidyverse) #packageVersion("tidyverse"); citation("tidyverse")
library(adegenet) #packageVersion("adegenet"); citation("adegenet")
library(WGCNA) #packageVersion("WGCNA"); citation("WGCNA")

countData <- read.table("B37_counts_newtranscriptome.txt")
head(countData)
length(countData[,1])
#26660

names(countData)=sub(".fastq.trim.sam.counts","",names(countData))
names(countData)
#row.names(countData)=sub("", "isogroup", rownames(countData))
names(countData)

#removing outlier samples (generally ones that failed to sequence)
countData$b37_10_7.5_4 <- NULL
countData$b37_24_7.5_2 <- NULL
countData$b37_24_7.5_4 <- NULL
countData$b37_4_8_2 <- NULL
countData$b37_4_8_4 <- NULL

totalCounts=colSums(countData)
totalCounts
barplot(totalCounts)
#    b37_0_8_1    b37_0_8_2    b37_0_8_3    b37_0_8_4 b37_10_7.5_1 
#       402011       437684       407762       126365       377236 
# b37_10_7.5_2 b37_10_7.5_3   b37_10_8_1   b37_10_8_2   b37_10_8_3 
#       522103       621347       223739       225729       260243 
#   b37_10_8_4 b37_24_7.5_1 b37_24_7.5_3   b37_24_8_1   b37_24_8_2 
#       587962       306236       360527       178560       554924 
#   b37_24_8_3   b37_24_8_4 b37_48_7.5_1 b37_48_7.5_2 b37_48_7.5_3 
#       139452       891746       609382      1235041       178275 
# b37_48_7.5_4   b37_48_8_1   b37_48_8_2   b37_48_8_3   b37_48_8_4 
#       297157       443670       917558       289920       620377 
#  b37_4_7.5_1  b37_4_7.5_2  b37_4_7.5_3  b37_4_7.5_4    b37_4_8_1 
#       605458       916263       314277       324889       926678 
#    b37_4_8_3 
#        94004 

min(totalCounts) #94,004
max(totalCounts)  # 1,235,041

#Remove day 0 data since it can't be compared between treatments (only pH 8 treatment available)
countData$b37_0_8_1 <- NULL
countData$b37_0_8_2 <- NULL
countData$b37_0_8_3 <- NULL
countData$b37_0_8_4 <- NULL

treat=c( "pH7.5", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8", "pH8", "pH8","pH7.5", "pH7.5", "pH7.5", "pH7.5","pH8", "pH8")
hour=as.factor(c( "10", "10", "10", "10", "10", "10", "10", "24", "24", "24", "24", "24", "24","48", "48", "48","48", "48", "48", "48", "48", "4","4", "4", "4", "4", "4"))

colData=data.frame(treat, hour)
colData4 = colData[colData$hour == 4,]
colData10 = colData[colData$hour == 10,]
colData24 = colData[colData$hour == 24,]
colData48 = colData[colData$hour == 48,]


dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat+hour+treat*hour ) 
head(dds)
dds4 = dds[,colData$hour == 4]
head(dds4)
dds10 = dds[,colData$hour == 10]
head(dds10)
dds24 = dds[,colData$hour == 24]
head(dds24)
dds48 = dds[,colData$hour == 48]
head(dds48)


design(dds4) = ~treat
design(dds10) = ~treat
design(dds24) = ~treat
design(dds48) = ~treat

dds<-DESeq(dds)
dds4<-DESeq(dds4)
dds10<-DESeq(dds10)
dds24<-DESeq(dds24)
dds48<-DESeq(dds48)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
head(dds)
res<- results(dds)
head(dds4)
res4<- results(dds4)
head(dds10)
res10<- results(dds10)
head(dds24)
res24<- results(dds24)
head(dds48)
res48<- results(dds48)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot")
plotDispEsts(dds4, main="Dispersion plot (Hour 4)")
plotDispEsts(dds10, main="Dispersion plot (Hour 10)")
plotDispEsts(dds24, main="Dispersion plot (Hour 24)")
plotDispEsts(dds48, main="Dispersion plot (Hour 48)")

###Assess Pco2 treatments: pH8 vs pH7.5
#Hour 4
respH75.4 <- results(dds4, contrast=c("treat","pH7.5","pH8"))
table(respH75.4$padj<.1)
table(respH75.4$padj<.05)
table(respH75.4$padj<.01)
# 0.1=2/25789
# 0.05=0/25789
# 0.01=0/25789
summary(respH75.4)
#LFC>0: 1, 0.0039%
#LFC<0: 1, 0.0039%
#Outliers: 0, 0%
#Low counts: 0, 0%

nrow(respH75.4[respH75.4$padj<0.05 & !is.na(respH75.4$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0

#scatterplot of logarithmic fold changes versus the mean of normalized counts
plotMA(respH75.4, main="pH8 vs pH7.5 (Hour 4)", ylim=c(-1,1))

results75.4 <- as.data.frame(respH75.4)
head(results75.4)

nrow(respH75.4[respH75.4$padj<0.1 & respH75.4$log2FoldChange > 0 & !is.na(respH75.4$padj),])
nrow(respH75.4[respH75.4$padj<0.1 & respH75.4$log2FoldChange < -0 & !is.na(respH75.4$padj),])
#UP in 7.5 (Hour 4): 1
#DOWN in 7.5 (Hour 4): 1

write.table(respH75.4, file="b37_7.5.4_2017.txt", quote=F, sep="\t")

cp4 <- read.table("b37_7.5.4_2017.txt")
head(cp4)

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

MAttPlot(cp4)
#orange dots = upregulated, blue dots = downregulated

#Hour 10
respH75.10 <- results(dds10, contrast=c("treat","pH7.5","pH8"))
table(respH75.10$padj<.1)
table(respH75.10$padj<.05)
table(respH75.10$padj<.01)
# 0.1=0/25814
# 0.05=0/25814
# 0.01=0/25814
summary(respH75.10)
#LFC>0: 0, 0%
#LFC<0: 0, 0%
#Outliers: 3, 0.012%
#Low counts: 0, 0%

nrow(respH75.10[respH75.10$padj<0.05 & !is.na(respH75.10$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0

plotMA(respH75.10, main="pH8 vs pH7.5 (Hour 10)", ylim=c(-1,1))

results75.10 <- as.data.frame(respH75.10)
head(results75.10)

nrow(respH75.10[respH75.10$padj<0.1 & respH75.10$log2FoldChange > 0 & !is.na(respH75.10$padj),])
nrow(respH75.10[respH75.10$padj<0.1 & respH75.10$log2FoldChange < -0 & !is.na(respH75.10$padj),])
#UP in 7.5 (Hour 10): 0
#DOWN in 7.5 (Hour 10): 0

write.table(respH75.10, file="b37_7.5.10_2017.txt", quote=F, sep="\t")

cp10 <- read.table("b37_7.5.10_2017.txt")
head(cp10)

MAttPlot(cp10)
#orange dots = upregulated, blue dots = downregulated

#Hour 24
respH75.24 <- results(dds24, contrast=c("treat","pH7.5","pH8"))
table(respH75.24$padj<.1)
table(respH75.24$padj<.05)
table(respH75.24$padj<.01)
# 0.1=0/25690
# 0.05=0/25690
# 0.01=0/25690
summary(respH75.24)
#LFC>0: 0, 0%
#LFC<0: 0, 0%
#Outliers: 1, 0.0039%
#Low counts: 0, 0%

nrow(respH75.24[respH75.24$padj<0.05 & !is.na(respH75.24$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #0

plotMA(respH75.24, main="pH8 vs pH7.5 (Hour 24)", ylim=c(-1,1))

results75.24 <- as.data.frame(respH75.24)
head(results75.24)

nrow(respH75.24[respH75.24$padj<0.1 & respH75.24$log2FoldChange > 0 & !is.na(respH75.24$padj),])
nrow(respH75.24[respH75.24$padj<0.1 & respH75.24$log2FoldChange < -0 & !is.na(respH75.24$padj),])
#UP in 7.5 (Hour 24): 0
#DOWN in 7.5 (Hour 24): 0

write.table(respH75.24, file="b37_7.5.24_2017.txt", quote=F, sep="\t")

cp24 <- read.table("b37_7.5.24_2017.txt")
head(cp24)

MAttPlot(cp24)
#orange dots = upregulated, blue dots = downregulated

#Hour 48
respH75.48 <- results(dds48, contrast=c("treat","pH7.5","pH8"))
table(respH75.48$padj<.1)
table(respH75.48$padj<.05)
table(respH75.48$padj<.01)
# 0.1=33/15575
# 0.05=22/15575
# 0.01=18/15575
summary(respH75.48)
#LFC>0: 18, 0.069%
#LFC<0: 15, 0.057%
#Outliers: 7, 0.027%
#Low counts: 10687, 41%

nrow(respH75.48[respH75.48$padj<0.05 & !is.na(respH75.48$padj),])  
# Num significantly differentially expressed genes excluding the no/low count genes   #22

#scatterplot of logarithmic fold changes versus the mean of normalized counts
plotMA(respH75.48, main="pH8 vs pH7.5 (Hour 48)", ylim=c(-1,1))

results75.48 <- as.data.frame(respH75.48)
head(results75.48)

nrow(respH75.48[respH75.48$padj<0.1 & respH75.48$log2FoldChange > 0 & !is.na(respH75.48$padj),])
nrow(respH75.48[respH75.48$padj<0.1 & respH75.48$log2FoldChange < -0 & !is.na(respH75.48$padj),])
#UP in 7.5 (Hour 48): 18
#DOWN in 7.5 (Hour 48): 15

write.table(respH75.48, file="b37_7.5.48_2017.txt", quote=F, sep="\t")

cp48 <- read.table("b37_7.5.48_2017.txt")
head(cp48)

MAttPlot(cp48)
#orange dots = upregulated, blue dots = downregulated

ddslogged = DESeq2::rlog(dds, blind = TRUE)
dds04logged = DESeq2::rlog(dds4, blind = TRUE)
dds10logged = DESeq2::rlog(dds10, blind = TRUE)
dds24logged = DESeq2::rlog(dds24, blind = TRUE)
dds48logged = DESeq2::rlog(dds48, blind = TRUE)

write.csv(assay(dds04logged), "dds04h_rlogged.csv")
write.csv(assay(dds10logged), "dds10h_rlogged.csv")
write.csv(assay(dds24logged), "dds24h_rlogged.csv")
write.csv(assay(dds48logged), "dds48h_rlogged.csv")
write.csv(assay(ddslogged), "ddslarvaeh_rlogged.csv")

##make the GO table for MWU
#Hour 4
cp4$isogroup=row.names(cp4)
go_input_7.5.4 = cp4 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.4)
colnames(go_input_7.5.4) <- c("gene", "pval")
head(go_input_7.5.4)
write.csv(go_input_7.5.4, file="b37_7.5.4_rec_GO.csv", quote=F, row.names=FALSE)

#Hour 10
cp10$isogroup=row.names(cp10)
go_input_7.5.10 = cp10 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.10)
colnames(go_input_7.5.10) <- c("gene", "pval")
head(go_input_7.5.10)
write.csv(go_input_7.5.10, file="b37_7.5.10_rec_GO.csv", quote=F, row.names=FALSE)

#Hour 24
cp24$isogroup=row.names(cp24)
go_input_7.5.24 = cp24 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.24)
colnames(go_input_7.5.24) <- c("gene", "pval")
head(go_input_7.5.24)
write.csv(go_input_7.5.24, file="b37_7.5.24_rec_GO.csv", quote=F, row.names=FALSE)

#Hour 48
cp48$isogroup=row.names(cp48)
go_input_7.5.48 = cp48 %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(isogroup, mutated_p_updown)
head(go_input_7.5.48)
colnames(go_input_7.5.48) <- c("gene", "pval")
head(go_input_7.5.48)
write.csv(go_input_7.5.48, file="b37_7.5.48_rec_GO.csv", quote=F, row.names=FALSE)

#GO enrichment analysis
#Hour 4

#######MF
# Edit these to match your data file names: 
input="b37_7.5.4_rec_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
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
#No significant GO terms

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

#GO enrichment analysis
#Hour 10

#######MF
input="b37_7.5.10_rec_GO.csv" 
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
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#GO enrichment analysis
#Hour 10

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

#GO enrichment analysis
#Hour 24

#######MF
# Edit these to match your data file names: 
input="b37_7.5.24_rec_GO.csv"
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
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results

#######BP
#No significant GO terms

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
                  treeHeight=0.5,
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results

#GO enrichment analysis
#Hour 48

#######MF
# Edit these to match your data file names: 
input="b37_7.5.48_rec_GO.csv"
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
                  level1=0.05,
                  level2=0.01,
                  level3=0.001,
                  txtsize=1.2,
                  treeHeight=0.5
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
                  level1=0.05,
                  level2=0.01,
                  level3=0.001,
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
                  level1=0.05,
                  level2=0.01,
                  level3=0.001,
                  txtsize=1.2,
                  treeHeight=0.5
)

# text representation of results, with actual adjusted p-values
results


################################################################################
#PCAs
PCA_colors <- ifelse(colData$treat == 'pH7.5', 'darkgoldenrod2', 
                     ifelse(colData$treat == 'pH8','dodgerblue3', 'darkolivegreen3'))

pt_sym <- ifelse(colData$hour == '4', 15, 
                 ifelse(colData$hour == '10',17,
                        ifelse(colData$hour == '24',16, 18)))
larv_pcadata = DESeq2::plotPCA(ddslogged, intgroup ="treat", returnData = TRUE)
larv_percentVar = round(100* attr(larv_pcadata, "percentVar"))
larv_pca =prcomp(t(assay(ddslogged)), center = TRUE, scale. = FALSE)

PCA_larv =as.data.frame(larv_pca$x)%>%
  dplyr::select(PC1, PC2)


plot(PCA_larv,choices=c(1,2),display="sites",type="n")
points(PCA_larv,col=PCA_colors, pch=pt_sym)
#ordispider(ccaa,groups=treat,col=c("green", "blue"))
#ordispider(ccaa,groups=conds$hour,col=c("green", "blue", "red", "yellow"))
ordihull(PCA_larv,groups=colData$hour,draw="polygon",col="grey95",label=T)
#ordihull(PCA_larv,groups=colData$treat,draw="polygon",border=PCA_colors,label=T)
legend(x="topleft", c("pH 7.5", "pH 8.0"), fill= c('darkgoldenrod2','dodgerblue3'))

adonis(larv_pca$x ~ colData$treat + colData$hour, method='eu')
#############################################################################
#rlog

dds04rlogged = read.table("dds04h_rlogged.csv")
dds04rlogged = as.data.frame(dds04rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7"), ",")
dds04rlogged$`1` = NULL
rownames(dds04rlogged) = dds04rlogged$V1
colnames(dds04rlogged) = dds04rlogged[1,]
colnames(dds04rlogged)[1] = "Gene"

dds10rlogged = read.table("dds10h_rlogged.csv")
dds10rlogged = as.data.frame(dds10rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7", "8"), ",")
dds10rlogged$`1` = NULL
rownames(dds10rlogged) = dds10rlogged$V1
colnames(dds10rlogged) = dds10rlogged[1,]
colnames(dds10rlogged)[1] = "Gene"

dds24rlogged = read.table("dds24h_rlogged.csv")
dds24rlogged = as.data.frame(dds24rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7"), ",")
dds24rlogged$`1` = NULL
rownames(dds24rlogged) = dds24rlogged$V1
colnames(dds24rlogged) = dds24rlogged[1,]
colnames(dds24rlogged)[1] = "Gene"

dds48rlogged = read.table("dds48h_rlogged.csv")
dds48rlogged = as.data.frame(dds48rlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7", "8", "9"), ",")
dds48rlogged$`1` = NULL
rownames(dds48rlogged) = dds48rlogged$V1
colnames(dds48rlogged) = dds48rlogged[1,]
colnames(dds48rlogged)[1] = "Gene"

ddsrlogged = read.table("ddslarvaeh_rlogged.csv")
ddsrlogged = as.data.frame(ddsrlogged)%>%
  separate(V2,c("1","2","3","4","5","6","7", "8", "9", "10", "11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"), ",")
ddsrlogged$`1` = NULL
rownames(ddsrlogged) = ddsrlogged$V1
colnames(ddsrlogged) = ddsrlogged[1,]
colnames(ddsrlogged)[1] = "Gene"

res4h <- results(dds4, contrast=c("treat","pH7.5","pH8"))
res4h = as.data.frame(res4h)
res4h = merge(res4h, dds04rlogged, by=0)
rownames(res4h) <- res4h$Gene

res10h <- results(dds10, contrast=c("treat","pH7.5","pH8"))
res10h = as.data.frame(res10h)
res10h = merge(res10h, dds10rlogged, by=0)
rownames(res10h) <- res10h$Gene

res24h <- results(dds24, contrast=c("treat","pH7.5","pH8"))
res24h = as.data.frame(res24h)
res24h = merge(res24h, dds24rlogged, by=0)
rownames(res24h) <- res24h$Gene

res48h <- results(dds48, contrast=c("treat","pH7.5","pH8"))
res48h = as.data.frame(res48h)
res48h = merge(res48h, dds48rlogged, by=0)
rownames(res48h) <- res48h$Gene

resh <- results(dds, contrast=c("treat","pH7.5","pH8"))
resh = as.data.frame(resh)
resh = merge(resh, ddsrlogged, by=0)
rownames(resh) <- resh$Gene

write.csv(res4h, "Hour4_Results.csv", row.names = TRUE)
write.csv(res10h, "Hour10_Results.csv", row.names = TRUE)
write.csv(res24h, "Hour24_Results.csv", row.names = TRUE)
write.csv(res48h, "Hour48_Results.csv", row.names = TRUE)
write.csv(resh, "48hLarvae_Results.csv", row.names = TRUE)
###################################################################
#DEG heatmaps 
library(tidyverse)
iso2gene = read.delim("crepidula_iso2gene.tab", header=FALSE)
rownames(iso2gene) <- iso2gene$V1


#Hour 48
Hour48_Results = read.csv("Hour48_Results.csv")
rownames(Hour48_Results) = Hour48_Results$X
Hour48_Results$X = NULL
rld_Hour48_data =as.data.frame(Hour48_Results[,c((6:7),(9:16))])

p.val=0.10 # FDR cutoff
conds_Hour48=rld_Hour48_data[rld_Hour48_data$padj<=p.val & !is.na(rld_Hour48_data$padj),]
length(conds_Hour48[,1])
#33
exp_Hour48=conds_Hour48[,3:10] #change to your columns
head(exp_Hour48)
means_Hour48=apply(exp_Hour48,1,mean) # means of rows
explc_Hour48=exp_Hour48-means_Hour48 # subtracting them


Hour48_df_all_iso <- explc_Hour48 %>%
  rownames_to_column("V1") %>%
  left_join(iso2gene) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))

Hour48_unanno=Hour48_df_all_iso[,2:10]

Hour48_df_only_anno <- Hour48_df_all_iso %>%
  filter(!is.na(V2))
rownames(Hour48_df_only_anno) <- make.unique(Hour48_df_only_anno$V2)
Hour48_df_only_anno$V1 = NULL

Hour48_df_unanno <- Hour48_df_all_iso %>%
  filter(is.na(V2))
rownames(Hour48_df_unanno) <- Hour48_df_unanno$V1
Hour48_df_unanno$V1 = NULL
Hour48_df_unanno$V2 = NULL

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(Hour48_df_only_anno)

Hour48_my_sample_col <- data.frame(Source = c("pH 7.5", "pH 7.5","pH 7.5", "pH 7.5", "pH 8", "pH 8", "pH 8","pH 8"))
Hour48_df_only_anno$V2 = NULL
row.names(Hour48_my_sample_col) = colnames(Hour48_df_only_anno)


my_colour_48h = list(
  Source = c('pH 7.5' = 'darkgoldenrod2', 'pH 8' = 'dodgerblue3'))

# big heat map of all annotated genes
pheatmap(Hour48_df_only_anno,cluster_cols=T,scale="row", color=col0, annotation_col= Hour48_my_sample_col, annotation_colors =my_colour_48h, show_rownames = T, show_colnames = F, border_color = "NA" )
dev.off()
# big heat map of all DEGs, no annotation too
pheatmap(Hour48_df_unanno,cluster_cols=T,scale="row",color=col0, annotation_col= Hour48_my_sample_col, annotation_colors =my_colour_48h, show_rownames = F,show_colnames = F, border_color = "NA" )
dev.off()

