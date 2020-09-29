#load libraries
library(devtools)
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
library(org.Mm.eg.db)
library(ggrepel)
library(knitr)
library(GOstats)
library(gage)
library(gageData)
library(Select)
library(kableExtra)
library(readxl)
library(ComplexHeatmap)

#create count table
counts <- read_excel("~/Documents/DrKwan/HKUlabdata/hkucounts.xlsx")
counts <- as.data.frame(counts)
counts <- counts[!duplicated(counts[,c("Symbol")]),]
row.names(counts) <- counts$Symbol
counts <- subset(counts, select = -c(1))
head(counts)

#create DGEList object 
d0 <- DGEList(counts)

#calculate normalization factors
d0 <- calcNormFactors(d0)
d0

#filter low-expressed genes (significantly reduces the number of genes in this case)
#cutoff <- 1
#drop <- which(apply(cpm(d0),1,max)<cutoff)
#d <- d0[-drop,]
#dim(d)

#derive experiment information from the sample names
snames <- colnames(counts)
snames
group <- substr(snames, 1, nchar(snames) -2)

plotMDS(d, col = as.numeric(group))

#specify the model to be fitted
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)

#fit linear models in limma
fit <- lmFit(y, mm)
head(coef(fit))
contr <- makeContrasts(groupMg - groupCTL, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)

#perform bayes smoothing
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

