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

## Example Dataset

Mitigating SOX2-potentiated Immune Escape of Head and Neck Squamous Cell Carcinoma with a STING-inducing Nanosatellite Vaccine
Yee Sun Tan et. al.
4 Samples

**Expected heatmap**
![nihms968654f1](https://user-images.githubusercontent.com/52653630/103473320-373abf80-4dd2-11eb-8219-230d322d6e6b.jpg)

## Notes

**large sample size and multiple layers**

DESeq2 compares two conditions/treatments. Therefore, as many experiments have more than two conditions, DESeq2 has to be run more than once. Furthermore, this reduces the number of columns for the heatmap and makes it a clean and readable visualisation








