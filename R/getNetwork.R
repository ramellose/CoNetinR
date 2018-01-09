#' @title Build Network
#' @description Build a network using Pearson, Spearman, Kullback-Leibler or Bray-Curtis
#'
#' @details Takes care of all steps of network construction, by also performing rarefication, normalization, permutation and bootstrapping.
#' If ReNorm and permutandboot are both set to TRUE, this is equal to the ReBoot procedure implemented in CoNet.
#' The original CoNet implementation is available at: http://psbweb05.psb.ugent.be/conet/
#'
#' @param mat user-provided matrix, can be left empty to generate one
#' @param method network construction method, values: "pearson", "spearman", "kld" or "bray"
#' @param T.up upper threshold for scores
#' @param T.down lower threshold for scores
#' @param shuffle.samples shuffle each sample before network construction
#' @param norm normalize matrix before network construction
#' @param rarefy if TRUE, rarefy matrix to the minimum total read count, if larger than 1, rarefy matrix to the given total read count
#' @param stand.rows standardize rows by dividing each entry by its corresponding row sum, applied after normalization and rarefaction and before network construction
#' @param pval.cor compute p-values of correlations with cor.test (cannot be selected together with permut, renorm and/or permutandboot or a method that is not a correlation)
#' @param permut compute p-values on edges with a simple permutation test
#' @param renorm compute p-values with a permutation test, using renormalization
#' @param permutandboot compute p-values from both permutation and bootstrap distribution
#' @param iters number of iterations for the permutation test
#' @param bh multiple-test-correct using Benjamini-Hochberg
#' @param min.occ only keep rows with at least the given number of non-zero values (carried out before network construction)
#' @param keep.filtered sum all filtered rows and add the sum vector as additional row
#' @param plot plot score or, if permut is true, p-value distribution
#' @param report.full in addition to positive edge percentage, return other output variables
#' @param verbose print the number of positive and negative edges and details of p-value computation
#' @return postive edge percentage and more statistics depending on report.full
#'
#' @importFrom stats cor.test p.adjust cor
#' @importFrom gdata lowerTriangle
#' @importFrom vegan vegdist
#' @importFrom graphics hist
#'
#' @export
getNetwork<-function(mat = matrix(), method="spearman", T.up=0.2, T.down=-0.2, shuffle.samples=F, norm=FALSE, rarefy=0, stand.rows=FALSE, pval.cor=F, permut=F, renorm=F, permutandboot=F, iters=100, bh=F, min.occ=0, keep.filtered=F, plot=F, report.full=F, verbose=F){
  min.occ.default=1
  keep.indices = c()
  if(pval.cor == TRUE & (method == "kld" || method == "bray")){
    stop("P-value computation with cor.test is only possible for correlations!")
  }
  if(pval.cor == TRUE & (permut == TRUE || permutandboot == TRUE || renorm == TRUE)){
    stop("P-value computation with cor.test cannot be combined with permut, renorm or permutandboot!")
  }
  # remove rows consisting of zeros only
  if(min.occ > 0){
    min.occ.default = min.occ
    for(i in 1:nrow(mat)){
      occ = 0
      for(j in 1:ncol(mat)){
        if(mat[i,j] > 0){
          occ = occ + 1
        }
      }
      if(occ >= min.occ.default){
        keep.indices = c(keep.indices,i)
      }
    }
    if(length(keep.indices) == 0){
      stop("No rows left in matrix after rare taxon filtering!")
    }
    is.summed=FALSE
    summed.row=c()
    if(keep.filtered == TRUE){
      filtered.rows=mat[setdiff(c(1:nrow(mat)),keep.indices),]
      if(nrow(filtered.rows) > 0){
        is.summed=TRUE
        summed.row = filtered.rows[1,]
        if(nrow(filtered.rows) > 1){
          for(f.index in 2:nrow(filtered.rows)){
            summed.row=summed.row+filtered.rows[f.index,]
          }
        }
      }
    }
    mat = mat[keep.indices,]
    if(keep.filtered == TRUE & is.summed == TRUE){
      mat=rbind(mat,summed.row)
    }
  }
  # if requested, rarefy matrix
  if(rarefy == TRUE){
    mat = seqtime::rarefyFilter(mat,min = 0)
  }else if(rarefy > 1){
    mat = seqtime::rarefyFilter(mat, min=rarefy)
  }
  # if requested, normalize matrix
  if(norm == TRUE){
    mat=seqtime::normalize(mat)
  }
  # normalize row-wise
  if(stand.rows == TRUE){
    mat = t(seqtime::normalize(t(mat)))
  }
  pvals = matrix(nrow=nrow(mat),ncol=nrow(mat))
  # compute p-values
  if(permut == TRUE || pval.cor == TRUE){
    if(iters < 1){
      stop("iters should be at least 1.")
    }
    for(i in 1:nrow(mat)){
      for(j in 1:i){
        if(i != j){
          if(pval.cor == TRUE){
            pvals[i, j] = cor.test(mat[i,],mat[j,],method=method)$p.value
          }else{
            pvals[i, j] = getPval(mat, i, j, method=method, N.rand=iters, renorm=renorm, permutandboot=permutandboot, verbose=verbose)
          }
          pvals[j, i] = pvals[i, j]
        }
      }
    }
    # if requested, carry out Benjamini-Hochberg correction
    if(bh == TRUE){
      pvec = pvals[lower.tri(pvals)]
      pvec=p.adjust(pvec,method="BH")
      # only lower triangle is used later on
      lowerTriangle(pvals)=pvec
    }
  }
  if(method == "spearman"){
    scores=cor(t(mat),method="spearman")
  }else if(method == "pearson"){
    scores=cor(t(mat),method="pearson")
  }else if(method == "bray"){
    scores = vegdist(mat, method="bray")
    scores=as.matrix(scores)
  }else if(method == "kld"){
    scores = computeKld(mat)
  }else{
    stop("Choose either spearman, pearson, kld or bray as a method.")
  }
  # if p-values were calculated, set all correlations with p-value above 0.05 to 0, set Bray Curtis scores to 0.5 and KLD scores to a value between the lower and upper threshold
  if(permut == TRUE || pval.cor == TRUE){
    FILTER=pvals>0.05
    if(method == "bray"){
      # Bray-Curtis dissimilarity of 0.5 is filtered out
      scores[FILTER]=0.5
    }else if(method == "kld"){
      # KLD value between the thresholds is filtered out
      scores[FILTER]=T.down+(T.up-T.down)/2
    }else{
      # correlation of zero is filtered out
      scores[FILTER]=0
    }
  }

  # get lower triangle of score matrix
  scores=scores[lower.tri(scores)]
  storage=scores
  # discard missing values
  scores=scores[!is.na(scores)]
  if(method == "bray" || method == "kld"){
    neg.edges=scores[scores>T.up]
    pos.edges=scores[scores<(T.down)]
  }else{
    pos.edges=scores[scores>T.up]
    neg.edges=scores[scores<(T.down)]
  }
  if(plot == T){
    values = scores
    what = "Score"
    if(permut == TRUE){
      values = pvals
      what = "P-value"
    }
    hist(values, main=paste(what," distribution for method ",method,sep=""), xlab=paste(what,"s",sep=""), ylab="Frequency")
  }
  num.pos=length(pos.edges)
  num.neg=length(neg.edges)
  total=num.pos+num.neg
  one.percent=total/100
  pos.percent=num.pos/one.percent
  neg.percent=num.neg/one.percent
  if(verbose == T){
    print(paste(num.pos,"positive edges",sep=" "))
    print(paste(num.neg,"negative edges",sep=" "))
    print(paste(pos.percent,"positive percent",sep=" "))
    print(paste(neg.percent,"negative percent",sep=" "))
  }
  if(report.full == TRUE){
    # append unfiltered score list
    out=list(pos.percent, neg.percent, num.pos, num.neg, total, storage, pvals)
    names(out)=c("pospercent", "negpercent","posnumber","negnumber","total","scores","pvalues")
    return(out)
  }else{
    return(pos.percent)
  }
}
