The following code worked on an older version of R

```
library(biomaRt)
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
gene_ids <- annotationFile[,6]
foo <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
             values = gene_ids,
             mart = ensembl)
colnames(foo) <- c("GeneID", "Gene_Symbol")
```

Currently there are the following errors

```
> library(biomaRt)
> ensembl <- useMart("ensembl")
> datasets <- listDatasets(ensembl)
> ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
Ensembl site unresponsive, trying useast mirror
Ensembl site unresponsive, trying asia mirror
> gene_ids <- annotationFile[,6]
> foo <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
+              values = gene_ids,
+              mart = ensembl)
Error: failed to load external entity "http://www.ensembl.org/info/website/archives/index.html?redirect=no"
> colnames(foo) <- c("GeneID", "Gene_Symbol")
Error in colnames(foo) <- c("GeneID", "Gene_Symbol") : 
  object 'foo' not found
```

Install developmental version from github 
```
BiocManager::install('grimbough/biomaRt')
```

Reference
https://support.bioconductor.org/p/134524/



