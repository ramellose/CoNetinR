---
title: "Tour through CoNet in R"
author: "Lisa Rottjers"
date: "`r Sys.Date()`"
output: 
  md_document    
    toc: true
vignette: >
  %\VignetteIndexEntry{seqtime examples: Properties of time series from different models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  \usepackage[utf8]{inputenc}
---


```{r setup, include=FALSE}
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_seqtime_tour/")
```

We start by loading the CoNetinR, vegan, seqtime and gdata libraries. 

```{r, message=FALSE, warning=FALSE}
library(CoNetinR)
library(vegan)
library(gdata)
library(seqtime)
```

The next step is to generate an interaction matrix, and from it, a dataset. 

```{r}
N = 50
S = 40
A = generateA(N, "klemm", pep=10, c =0.05)
dataset = generateDataSet(S, A)
```

Next, we can use the CoNet adaptation to try and find back the original interaction matrix. We are going to use the Spearman method, and first get the Spearman scores. 

```{r}
scores = getNetwork(mat = A, method="spearman", T.up=0.2, T.down=-0.2, shuffle.samples=F, norm=TRUE, rarefy=0, stand.rows=F, pval.cor=F, permut=F, renorm=F, permutandboot=T, iters=100, bh=T, min.occ=0, keep.filtered=F, plot=F, report.full=T, verbose=F)
scores = scores$scores

```

We also need to get the p-values. 

```{r, fig.height = 6, fig.width = 6}
pmatrix = getNetwork(mat = dataset, method="spearman", T.up=0.2, T.down=-0.2, shuffle.samples=F, norm=TRUE, rarefy=0, stand.rows=F, pval.cor=T, permut=F, renorm=F, permutandboot=F, iters=100, bh=T, min.occ=0, keep.filtered=F, plot=F, report.full=T, verbose=F)
pmatrix = pmatrix$pvalues

```

Of course, now we have the Spearman correlations and the p-values. We can turn that into an adjacency matrix.

```{r}
adjmatrix = matrix(nrow = N, ncol = N)
adjmatrix[lower.tri(adjmatrix)] = scores
adjmatrix = t(adjmatrix)
adjmatrix[lower.tri(adjmatrix)] = scores
for (i in 1:N){
   for (j in 1:N){
    if (is.na(adjmatrix[i,j])){
       adjmatrix[i,j] = 0
    }
    else if (pmatrix[i,j] > 0.05){
       adjmatrix[i,j] = 0
    }
   }
 }
```

The adjacency matrix can be plotted with seqtime. The first figure is the original interaction matrix, while the second is the inferred matrix.

```{r}
plotA(A)
plotA(adjmatrix)
```

