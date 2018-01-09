#' @title Compute KLD
#' @description Compute Kullback-Leibler dissimilarity
#'
#' @details Outputs a Kullback-Leibler dissimilarity (symmetrized divergence) matrix.Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#'
#' @param x matrix-like object or vector with non-negative numbers
#' @param pseudocount pseudocount = this value is added to each zero
#
#' @return kullback-leibler dissimilarity matrix
#'
#' @importFrom stats cor.test p.adjust cor
#' @importFrom vegan vegdist
#' @importFrom graphics hist
#'
#' @export
computeKld=function(x, pseudocount=0.00000001){
  # diagonal is zero
  kld=matrix(data=0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:i){
      kld[i,j]=get.kld(x[i,],x[j,], pseudocount=pseudocount)
      kld[j,i]=kld[i,j]
    }
  }
  kld
}

get.kld=function(x,y, pseudocount=0.00000001){
  if(length(x) != length(y)){
    stop("The two vectors should have the same length!")
  }
  x[x==0]=pseudocount
  y[y==0]=pseudocount
  dis = 0
  x = x/sum(x)
  y = y/sum(y)
  for(i in 1:length(x)){
    if(!is.nan(x[i]) && !is.nan(y[i])){
      ratioxy = log(x[i]/y[i])
      ratioyx = log(y[i]/x[i])
      dis = x[i]*ratioxy+y[i]*ratioyx + dis
    }
  }
  dis
}
