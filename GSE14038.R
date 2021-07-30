---
title: "GSE14038"
author: "Teresa Rozza"
date: "7/22/2021"
output: html_document
---

```{r include=FALSE}
require(knitr)
opts_chunk$set(
concordance=FALSE, echo=TRUE,  warning=FALSE, error=FALSE, message=FALSE)
```

```{r}
# Instalar Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) # install.packages("BiocManager")
BiocManager::install(version = "3.13")
# Instalar y cargar libraries
if (!requireNamespace("GEOquery", quietly = TRUE))
BiocManager::install("GEOquery") 
if (!requireNamespace("limma", quietly = TRUE))
BiocManager::install("limma") 
if (!requireNamespace("umap", quietly = TRUE))
BiocManager::install("umap") 
if (!requireNamespace("pvca", quietly = TRUE))
BiocManager::install("pvca")
if (!requireNamespace("oligo", quietly = TRUE))
BiocManager::install("oligo")
if (!requireNamespace("pd.mogene.2.1.st", quietly = TRUE))
BiocManager::install("pd.mogene.2.1.st")
if (!requireNamespace("Biobase", quietly = TRUE))
BiocManager::install("Biobase")
if (!requireNamespace("genefilter", quietly = TRUE))
BiocManager::install("genefilter")
if (!requireNamespace("mogene21sttranscriptcluster.db", quietly = TRUE))
BiocManager::install("mogene21sttranscriptcluster.db")
if (!requireNamespace("annotate", quietly = TRUE))
BiocManager::install("annotate")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
BiocManager::install("org.Mm.eg.db")
if (!requireNamespace("ReactomePA", quietly = TRUE))
BiocManager::install("ReactomePA")
if (!requireNamespace("reactome.db", quietly = TRUE))
BiocManager::install("reactome.db")
```

```{r}
library(GEOquery)
library(limma)
library(umap)
library(pvca)
library(oligo)
library(pd.mogene.2.1.st)
library(Biobase)
library(genefilter)
library(mogene21sttranscriptcluster.db)
library(annotate)
library(org.Mm.eg.db)
library(ReactomePA)
library(reactome.db)
```


```{r}
FULL_DIR <- "/Users/teresa/Desktop/hack4rare/GSE14038_RAW/"
RAW_DIR <- "~/Desktop/hack4rare/GSE14038_RAW"
datadir <- setwd(RAW_DIR)
```


# STEP 1 - Obtain Microarray Data

```{r}
# Store the dataset ids in a vector GEO_DATASETS just in case you want to loop through several GEO ids
GEO_DATASETS <- "GSE14038"
```

```{r}
#loading dataset
gset <- getGEO("GSE14038", destdir=datadir, GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL7869", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) }
```


```{r}
#quality data
#library(arrayQualityMetrics)
#arrayQualityMetrics(gset)
```
From the file index.htlm we can see that only 4 specimen are marked once, so we can keep them in the analysis.


```{r}
#explore structure of data
gset
```

```{r}
head(exprs(gset))
```


```{r}
# Column-wise summary statistics,3
round(apply(ex,2, summary),3)  
```


```{r}
# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE14038")
```


```{r}
boxplot(ex)
```


```{r}
colnames(pData(gset))
```
```{r}
pData(gset)$data_processing[1]
```

We now have the metadata, phenodata (sample data), and experimental data (microarray intensities) for this experiment.Now we have to downloads the raw data files to your computer. For Affymetrix microarrays, the raw data format is a CEL file.

```{r}
#downloading raw data
celFiles <- list.celfiles(RAW_DIR, full.names = TRUE)
Data <- read.celfiles(celFiles)
```
```{r}
head(Data)
```

# STEP 2 - Processing Microarray Data

The first step to process CEL files with "affy" is to tell the program the directory that it can find the CEL files and the corresponding file with phenoData (sample data).

## Preparing the Phenodata

The rownames of the data frame must be the same as the CEL file names

```{r}
head(celFiles)
```

```{r}
#eliminate unseful text
celFiles <- gsub(FULL_DIR, "", celFiles)
```


```{r}
my.pdata <- as.data.frame(pData(gset), stringsAsFactors=F)
my.pdata <- my.pdata[, c("title", "geo_accession", "description")]
my.pdata <- my.pdata[order(rownames(my.pdata)), ]
head(my.pdata, 10)
```
```{r}
group <- my.pdata$description
```


```{r}
#eliminate unseful text
group <- gsub("Gene expression data from", "", group)
my.pdata <- cbind(my.pdata, group)
my.pdata
```


```{r}
head(celFiles)
```

```{r}
head(rownames(my.pdata))
```


```{r}
table(rownames(my.pdata) == celFiles)
```

```{r}
temp.rownames <- paste(rownames(my.pdata), ".CEL", sep="")
table(temp.rownames == celFiles)
```


```{r}
rownames(my.pdata) <- temp.rownames
rm(temp.rownames)
table(rownames(my.pdata) == celFiles)
```


```{r}
head(my.pdata)
```

## Reading the CEL Files

Now that we have directory of CEL files and a corresponding data frame with the phenoData, we can read the CEL files into R

```{r}
list.files(RAW_DIR)
```


```{r}
##Perform affy normalization
library(affy)
my.affy <- ReadAffy(celfile.path = RAW_DIR, phenoData = my.pdata)
show(my.affy)
```

```{r}
head(exprs(my.affy),5)
```

```{r}
dim(exprs(my.affy))
```

```{r}
colnames(pData(my.affy))
```


## Calculating Gene Expression Measurements

the expression data consists of 1354896 individual probe intensities for each of our 86 samples. We need to combine the individual probe intensities to probeset-level (gene-level) measurements. Each probeset typically consists of 11 - 20 individual probes. We will use the rma() function to combine the individual probe intensities to a probeset intensity. 

```{r}
##Calculate gene level expression measures
my.rma <- rma(my.affy, normalize=F, background=F)
head(exprs(my.rma))
```

```{r}
pData(my.rma)
```


```{r}
# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE14038", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)
```


## STEP 3 - EXPLORATORY ANALYSIS 
```{r}
# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
```


```{r pca}
#computing prncipal components and loadings.
pcX<-prcomp(t(ex), scale=TRUE) 
loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
```


```{r}
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, xlim=c(-150, 150))
title("Principal components (PCA)")
text(pcX$x[,1],pcX$x[,2], colnames(ex), pos=4)
```



```{r}
sds <- apply(ex, 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
```

## STEP 4 - Differential gene Expression

```{r}
plotMDS(exprs(my.rma), labels=pData(my.rma)$group, top=500, gene.selection="common", main="MDS Plot to Compare Replicates")
```

Alternatively a hierarchichal clustering can be applied to detect any expected (or unexpected grouping of the samples).


```{r}
cluster.dat <- exprs(my.rma)
gene.mean <- apply(cluster.dat, 1, mean)
gene.sd <- apply(cluster.dat, 1, sd)
cluster.dat <- sweep(cluster.dat, 1, gene.mean, "-")
cluster.dat <- sweep(cluster.dat, 1, gene.sd, "/")
my.dist <- dist(t(cluster.dat), method="euclidean")
my.hclust <- hclust(my.dist, method="average")
my.hclust$labels <- pData(my.rma)$group
plot(my.hclust, cex=0.75, main="Comparison of Biological Replicates", xlab="Euclidean Distance")
```


## Design of the experiment


```{r}
##determine the average effect (coefficient) for each condition
library(limma)
designMat<- model.matrix(~0+group, pData(my.rma))
colnames(designMat) <- c("neuro_whole", "neuro_cell", "schwann_cell", "mouse_whole")
print(designMat)
```

## Fitting coefficients

"neuro_whole", "neuro_cell", "schwann_cell", "mouse_whole"

```{r}
##specify the contrast of interest using the levels from the design matrix
cont.matrix <- makeContrasts (neuroVSschawnn_cell = neuro_cell-schwann_cell,
                              neuroVSmouse_whole = neuro_whole-mouse_whole,
                              neurowVSschwann = neuro_whole-schwann_cell,
                              INT = (neuro_cell-schwann_cell) - (neuro_whole-schwann_cell),
                              levels=designMat)
print(cont.matrix)
```



Now that we have a design matrix, we need to estimate the coefficients. For this design, we will essentially average the replicate arrays for each sample level. In addition, we will calculate standard deviations for each gene, and the average intensity for the genes across all microarrays.

```{r}
##determine the average effect (coefficient) for each treatment
my.fit <- lmFit(my.rma, designMat)
my.fit
```

## Obtaining lists of differentially expressed genes and performing statistical analysis

The eBayes function performs the tests, and there are several parameters (arguments) that can be changed. We are going to change only the proportion of genes that we expect to be differentially expressed. 
```{r}
#linear model fit
fit.main<-contrasts.fit(my.fit, cont.matrix)
names(fit.main)
```


```{r}
fit.Bayes<-eBayes(fit.main, proportion=0.1, trend=FALSE, robust=FALSE)
```


The `limma` package implements function `topTable` which contains, for a given contrast a list of genes ordered from smallest to biggest p--value which can be considered to be most to least differential expressed. For each gene the following statistics are provided:

- `logFC`: Mean difference between groups.  
-  `AveExpr`: Average expression of all genes in the comparison.
-  `t` : Moderated t-statistic (t-test-like statistic for the comparison).
-  `P.Value`: Test p--value.  
-  `adj.P.Val`: Adjusted p--value      
-  `B`: B-statistic: Posterior log odds of the gene of being vs non being differential expressed.
topTable will adjust the p-values and return the top genes that meet the cutoffs.

 ```{r}
#first comparison
topTab_neuroVSschawnn_cell <- topTable (fit.Bayes, number=nrow(fit.Bayes), coef="neuroVSschawnn_cell", adjust="fdr") 
head(topTab_neuroVSschawnn_cell)
```

```{r}
#second comparison
topTab_neurowVSschwann <- topTable (fit.Bayes, number=nrow(fit.Bayes), coef="neurowVSschwann", adjust="fdr") 
head(topTab_neurowVSschwann)
```

```{r}
#third comparison
topTab_neuroVSmouse_whole <- topTable (fit.Bayes, number=nrow(fit.Bayes), coef="neuroVSmouse_whole", adjust="fdr") 
head(topTab_neuroVSmouse_whole)
```

```{r}
#third comparison: Genes that behave differently between comparison 1 and 2
topTab_INT <- topTable(fit.Bayes, number=nrow(fit.Bayes), coef="INT", adjust="fdr") 
head(topTab_INT)
```


Finally, decideTests will make calls for DEGs by adjusting the p-values and applying a logFC cutoff similar to topTable.

```{r}
contrast.tests <- decideTests(fit.Bayes, method="separate", adjust.method="BH", p.value=0.05, lfc=0)
head(contrast.tests)
```
                                
                                
                                
## Assessing the Results
There are several ways to visualize the statistical results from the DGE analysis. Limma has a volcanoplot function.


```{r}
volcanoplot(fit.Bayes, coef=1, highlight=4, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
  abline(v=c(-1,1))
```

```{r}
volcanoplot(fit.Bayes, coef=2, highlight=4, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[2], sep="\n"))
  abline(v=c(-1,1))
```

```{r}
volcanoplot(fit.Bayes, coef=3, highlight=4, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[3], sep="\n"))
  abline(v=c(-1,1))
```

```{r}
volcanoplot(fit.Bayes, coef=4, highlight=4, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[4], sep="\n"))
  abline(v=c(-1,1))
```
                                
## Biological significance

```{r}
library (hgu133plus2.db)
columns(hgu133plus2.db)
```

```{r}
gene.data1 <- select(hgu133plus2.db, keys=rownames(topTab_neuroVSschawnn_cell), keytype="PROBEID", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
head(gene.data1)
```

```{r}
gene.data2 <- select(hgu133plus2.db, keys=rownames(topTab_neurowVSschwann), keytype="PROBEID", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
head(gene.data2)
```

```{r}
gene.data3 <- select(hgu133plus2.db, keys=rownames(topTab_neuroVSmouse_whole), keytype="PROBEID", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
head(gene.data3)
```


```{r}
gene.data4 <- select(hgu133plus2.db, keys=rownames(topTab_INT), keytype="PROBEID", columns=c("ENTREZID", "GENENAME", "SYMBOL"))
head(gene.data4)
```

                                
                                


