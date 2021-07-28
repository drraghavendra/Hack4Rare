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
BiocManager::install("GEOquery") 
BiocManager::install("limma") 
BiocManager::install("umap") 
BiocManager::install("pvca")
BiocManager::install("oligo")
BiocManager::install("pd.mogene.2.1.st")
BiocManager::install("Biobase")
BiocManager::install("genefilter")
BiocManager::install("mogene21sttranscriptcluster.db")
BiocManager::install("annotate")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("ReactomePA")
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
datadir <- setwd("~/Desktop/hack4rare/GSE14038_RAW")
```


## STEP 1 - LOADING AND NORMALIZING DATASET
```{r}
#loading dataset
gset <- getGEO("GSE14038", GSEMatrix =TRUE, getGPL=TRUE)
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
library(arrayQualityMetrics)
arrayQualityMetrics(gset)
```
From the file index.htlm we can see that only 4 specimen are marked once, so we can keep them in the analysis. 


```{r}
# Column-wise summary statistics,3
round(apply(ex,2, summary),3)  
```

```{r}
boxplot(ex)
```


```{r}
# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE14038", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)
```

```{r}
# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE14038")
```

## STEP 2 - EXPLORATORY ANALYSIS 
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

Alternatively a hierarchichal clustering can be applied to detect any expected (or unexpected grouping of the samples).


```{r codedendrogramcomputeHC}
clust.euclid.average <- hclust(dist(t(ex)),method="ward.D2")
```

```{r plotdendrograms, fig=T}
plot(clust.euclid.average, hang=-1)
```


```{r}
sds <- apply(ex, 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
```

## STEP 3 - DESIGN THE MODEL




