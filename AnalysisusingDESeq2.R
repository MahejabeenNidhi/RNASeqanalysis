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

#upload count matrix
count <- read.csv("~/Documents/DrKwan/Workshop/DESeq.csv")

#create data frame 
countData <- as.data.frame(countdata)
countData <- countData[!duplicated(countData[,c("Symbol")]),]
row.names(countData) <- countData$Symbol
countData <- countData[,-c(1:6)]

#data exploration
dim(countData)

columns = c("CTL_1", "CTL_2", "CTL_3", "Mg_3")
colnames(countData) = columns

hist(countData[,1], br=200, xlab="Number of Reads Counts per Feature", main="Histogram of Read Counts")

logCountData = log2(1+countData)
par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns
hist(logCountData[,1], main="Histogram of Log Read Counts", xlab="Log transformed counts")
boxplot(logCountData,las=3, main="Boxplot of Log Read Counts")

x <- logCountData
myColors = rainbow(dim(x)[2])
plot(density(x[,1]),col = myColors[1], lwd=2,
     xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
     ylim=c(0, max(density(x[,1])$y)+0.5) )

for( i in 2:dim(x)[2] )
  lines(density(x[,i]),col=myColors[i], lwd=2)
legend("topright", cex=0.5,colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	

detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}

groups = as.character ( detectGroups( colnames( countData ) ) )
groups

cell = c("CTL", "CTL", "CTL", "Mg")

colData = cbind(colnames(countData),cell)
colData

str(colData)

dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design= ~ cell)   # note that the study design is changed.
dds = DESeq(dds)  # main function
nrow(dds)

dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)

rld <- rlog(dds, blind = FALSE)

vsd <- vst(dds, blind = FALSE)

dds <- estimateSizeFactors(dds)

slog <- log2(counts(dds, normalized=TRUE)+1)

par(mfrow = c(1, 3))  # 3 columns
plot(slog[,1],slog[,2])
plot(assay(rld)[,1],assay(rld)[,2])
plot(assay(vsd)[,1],assay(vsd)[,2])

par(mfrow = c(1, 3))  # 3 columns
slog <- log2(counts(dds, normalized=TRUE)+1)
plot(slog[,1],slog[,2])
slog <- log2(counts(dds, normalized=TRUE)+4)
plot(slog[,1],slog[,2], xlim=c(0,20))
slog <- log2(counts(dds, normalized=TRUE)+20)
plot(slog[,1],slog[,2], xlim=c(0,20))

df <- bind_rows(
  as_data_frame(slog[,1:2]) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

plotPCA(rld, intgroup = c("cell")) + theme(aspect.ratio=1)

pca.object <- prcomp(t(assay(rld))) # PCA 
pcaData = as.data.frame(pca.object$x[,1:2]); 
pcaData = cbind(pcaData,detectGroups(colnames(assay(rld)) ))
colnames(pcaData) = c("PC1", "PC2", "Type")
percentVar=round(100*summary(pca.object)$importance[2,1:2],0)
#plot
p=ggplot(pcaData, aes(PC1, PC2, color=Type, shape = Type)) + geom_point(size=5) 
p=p+xlab(paste0("PC1: ",percentVar[1],"% variance")) 
p=p+ylab(paste0("PC2: ",percentVar[2],"% variance")) 
p=p+ggtitle("Principal component analysis (PCA)")+coord_fixed(ratio=1.0)+ 
  theme(plot.title = element_text(size = 16,hjust = 0.5)) + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text( size = 16),
        axis.text.y = element_text( size = 16),
        axis.title.x = element_text( size = 16),
        axis.title.y = element_text( size = 16) ) +
  theme(legend.text=element_text(size=16))
print(p)

dist2 <- function(x, ...)   # distance function = 1-PCC (Pearson's correlation coefficient)
  as.dist(1-cor(t(x), method="pearson"))

fit = cmdscale( dist2(t(assay(rld))) , eig=T, k=2)
mdsData <- as.data.frame(fit$points[,1:2]); 
mdsData <- cbind(mdsData,detectGroups(colnames(assay(rld))) )
colnames(mdsData) = c("x1", "x2", "Type")

p<-ggplot(mdsData, aes(x1, x2, color=Type, shape = Type)) + geom_point(size=5) 
p=p+xlab("Dimension 1") 
p=p+ylab("Dimension 2") 
p=p+ggtitle("Multidimensional scaling (MDS)")+ coord_fixed(ratio=1.)+ 
  theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
  theme(axis.text.x = element_text( size = 15),
        axis.text.y = element_text( size = 15),
        axis.title.x = element_text( size = 15),
        axis.title.y = element_text( size = 15) ) +
  theme(legend.text=element_text(size=15))
print(p)

library(gplots)

hclust2 <- function(x, method="average", ...)  # average linkage in hierarchical clustering
  hclust(x, method=method, ...)

n=100 # number of top genes by standard deviation

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

groups = detectGroups(colnames(x) )
groups.colors = rainbow(length(unique(groups) ) )


#pretty heat map
pheatmap(x)


