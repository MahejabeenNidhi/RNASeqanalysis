#load libraries
library(edgeR)
library(limma)
library(DESeq2)
library(Biobase)
library(EnhancedVolcano)
library(GenomicAlignments)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(viridisLite)
library(genefilter)
library(ggplot2)
library(readr)
library(dplyr)
library(gplots)
library(org.Hs.eg.db)
library(ggrepel)
library(knitr)
library(GOstats)
library(gage)
library(gageData)
library(Select)
library(kableExtra)
library(readxl)
library(ComplexHeatmap)

#load count matrix 
count <- read.csv("~/Documents/DrKwan/Workshop/DESeqwTrimmedData.csv")

#create data frame 
countData <- as.data.frame(count)
countData <- countData[!duplicated(countData[,c("Gene_Symbol")]),]
row.names(countData) <- countData$Gene_Symbol
countData <- countData[,c(7:17)]

#create column data

colData <- DataFrame(
  cell = c("SPtetramerCD8TCell", "SPtetramerCD8TCell", "SPtetramerCD8TCell", "SPtetramerCD8TCell", "SPtetramerCD8TCell", "Ly49NCD8TCell", "Ly49PCD8TCell", "Ly49NCD8TCell", "Ly49PCD8TCell", "Ly49NCD8TCell", "Ly49PCD8TCell"),
  treatment = c("MOGSP", "MOGSP", "MOG", "MOG", "MOG", "MOGSP", "MOGSP", "MOGSP", "MOGSP", "MOGSP", "MOGSP"),
  row.names = colnames(countData)
)

colData

#DESeq2

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design= ~ treatment + cell)   

dds$cell <- dropEmptyLevels(dds$cell)
dds$treatment <- dropEmptyLevels(dds$treatment)

#remove genes with low counts 

keepCounts <- rowSums(counts(dds)) >= 250
dds <- dds[keepCounts,]

#change the levels of "cell" 

levels(dds$cell)
levels(dds$cell) <- sub("-.*", "", levels(dds$cell))
levels(dds$cell)

#DESeq2 
#Interactions

#If the comparisons of interest are, for example, the effect of a condition for different sets of samples, a simpler approach than adding interaction terms explicitly to the design formula is to perform the following steps:
##combine the factors of interest into a single factor with all combinations of the original factors
##change the design to include just this factor, e.g. ~ group

dds$group <- factor(paste0(dds$treatment, dds$cell))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast = c("group", "MOGSPLy49PCD8TCell", "MOGSPLy49NCD8TCell"))

#p-values and adjusted p-values

resOrdered <- res[order(res$pvalue),]
summary(res)

#Benjamini-Hochberg two tailed adjusted P <0.005 from Wald's test
#Log2 Fold Change >0.75

res005 <- results(dds, alpha=0.0025, lfcThreshold = 0.75)
summary(res005)
sum(res005$padj < 0.0025, na.rm=TRUE)


#If the variable is continuous or an interaction term
#then the results can be extract name is one of elements returned by `resultsNames(dds)`.

#transform count data

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#Wald test individual steps

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#p-values and adjusted p-values

resOrderedCellTreatment <- res[order(res$pvalue),]
summary(resOrderedCellTreatment)

#more information on results columns

mcols(resOrderedCellTreatment)$description

df <- bind_rows(
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

dist2 <- function(x, ...)   # distance function = 1-PCC (Pearson's correlation coefficient)
  as.dist(1-cor(t(x), method="pearson"))

library(gplots)

hclust2 <- function(x, method="average", ...)  # average linkage in hierarchical clustering
  hclust(x, method=method, ...)

n=130 # number of top genes by standard deviation

x = assay(rld)
if(n>dim(x)[1]) n = dim(x)[1] # max	as data

x = x[order(apply(x,1,sd),decreasing=TRUE),]  # sort genes by standard deviation

x = x[1:n,]   # only keep the n genes

# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
x[x< cutoff] <- cutoff

#pretty heat map
pheatmap(x, fontsize_row = 3)
