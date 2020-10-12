## RNASeqanalysis

Create visualization to identify the genes that have been upregulated or downregulated and Gene Set Enrichment Analysis

## Dependency

To run this code the following R/Bioconductor packages are required:
* edgeR
* limma
* DESeq2
* Biobase
* gplots
* ggrepel
* gage
* Select
* ggplot2
* tidyverse
* genefilter
* EnhancedVolcano
* pheatmap
* org.Hs.eg.db
* org.Mm.eg.db

**Generating GeneCounts**

If faced with hardware limitations, run RNA STAR on Galaxy and select GeneCounts as part of the output. 

## Using the code

**Processing GeneCounts for DESeq2**

**Running DESeq2**

**Heatmap**

**EnhancedVolcano**

**Running limma for HTSanalyzeR2**

**HTSanalyzeR2**

## Notes

**large sample size and multiple layers**

DESeq2 compares two conditions/treatments. Therefore, as many experiments have more than two conditions, DESeq2 has to be run more than once. Furthermore, this reduces the number of columns for the heatmap and makes it a clean and readable visualisation

In the case of the mouse example,

|Sample|Cell|Treatment|
|------|----|---------|
|1|SPtetramer|MOGSP|
|2|SPtetramer|MOGSP|
|3|SPtetramer|MOG|
|4|SPtetramer|MOG|
|5|SPtetramer|MOG|
|6|Ly49-|MOGSP|
|7|Ly49+|MOGSP|
|8|Ly49-|MOGSP|
|9|Ly49+|MOGSP|
|10|Ly49-|MOGSP|
|11|Ly49+|MOGSP|

DESeq2 could be run twice

```
dds = DESeqDataSetFromMatrix(counts,colData=colData,~Ly +treatment)
dds = DESeq(dds)
```

Ly49Minus vs Ly49Plus

```
head(results(dds,c("Ly","Ly49plus","Ly49minus")))
log2 fold change (MLE): Ly Ly49plus vs Ly49minus 
Wald test p-value: Ly Ly49plus vs Ly49minus 
```

MOG vs MOGSP

```
results(dds,c("treatment","MOGSP","MOG"))
log2 fold change (MLE): treatment MOGSP vs MOG 
Wald test p-value: treatment MOGSP vs MOG 
```

As there are only 11 samples, all samples and their biological replicates could be included in the heatmap. 

For the human example, there were 30 samples, so including the biological replicates would make the visualisation difficult to interpret. 
Therefore, the matrix to feed into the heatmap is the foldchanges for 2 timepoints times 3 conditions. 
extract log fold change for 
* 6h xRNA vs 6h noRNA
* 16h xRNA vs 16h noRNA
* 6h pUUC vs 6h noRNA
* 16h pUUC vs 16h noRNA
* 6h R848 vs 6h vehicle
* 16h R848 vs 16h vehicle





