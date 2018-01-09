#' @title Get p-value
#' @description Get permuation-based p-value for association between two vectors.
#'
#' @details # Compute the association between two vectors using the given method and
#' compute its p-value using a permutation test. This method was adapted from R code by Fah Sahtirapongsasuti.
#' This method was adapted from CCREPE: http://huttenhower.sph.harvard.edu/ccrepe.
#' Emma Schwager et al Detecting statistically significant associtations between sparse and high dimensional compositional data. (In progress).
#'
#' @param matrix input matrix
#' @param x.index index of first vector in the input matrix
#' @param y.index index of second vector in the input matrix
#' @param N.rand number of iterations used for the permutation test
#' @param method similarity measure (supported measures are: "pearson", "spearman", "bray" and "kld")
#' @param renorm renormalize after permutation
#' @param permutandboot compute a bootstrap distribution in addition to the permutation distribution and
#                 return the p-value as the mean of the permutation distribution under the bootstrap distribution
#' @param plot plot the histogram of the permutation and, if permutandboot is true, of the bootstrap distribution
#' @param verbose print distribution properties and p-value
#'
#' @return p-value of the association
#'
#' @importFrom stats cor na.omit pnorm sd
#' @importFrom grDevices rgb
#' @importFrom seqtime normalize
#' @importFrom vegan vegdist
#' @importFrom graphics hist abline legend
#'
#' @export
getPval = function(matrix, x.index, y.index, N.rand=1000, method="spearman", renorm=F, permutandboot=F, plot=F, verbose=F) {
  x = matrix[x.index,]
  y = matrix[y.index,]
  lower.tail = TRUE
  # bray and kld are dissimilarities, so one-sided p-value needs to be computed from the upper tail
  if(method == "bray" || method == "kld"){
    lower.tail = FALSE
  }
  if(method == "spearman"){
    this.sim = cor(x, y, use="complete.obs", method="spearman")
  }else if(method == "pearson"){
    this.sim = cor(x, y, use="complete.obs", method="pearson")
  }else if(method == "bray"){
    this.sim= vegdist(rbind(x,y),method="bray")
  }else if(method == "kld"){
    this.sim=get.kld(x,y)
  }else{
    stop("Select either spearman, pearson, kld or bray as method.")
  }
  rand.sim = rep(NA, N.rand)
  boot.sim = rep(NA, N.rand)
  for (i in 1:N.rand) {
    rand = sample(x, length(x))
    if(renorm == T){
      mat.copy=matrix
      mat.copy[x.index,]=rand
      mat.copy = normalize(mat.copy)
      rand = mat.copy[x.index,]
      y = mat.copy[y.index,]
    }
    if(method == "spearman"){
      rand.sim[i] = cor(rand, y, method="spearman", use="complete.obs")
    }else if(method == "pearson"){
      rand.sim[i] = cor(rand, y, method="pearson",use="complete.obs")
    }else if(method == "bray"){
      rand.sim[i] = vegdist(rbind(rand,y),method="bray")
    }else if(method == "kld"){
      rand.sim[i] = computeKld(rand,y)
    }
  }
  rand.sim = na.omit(rand.sim)
  if(plot == T){
    col1=rgb(0,0,1,1/3)
    col2=rgb(1,0,0,1/3)
    hist(rand.sim,col=col1)
    abline(v=mean(rand.sim),col="blue")
  }
  if(permutandboot){
    x=matrix[x.index,]
    y=matrix[y.index,]
    for (i in 1:N.rand) {
      rand.idx = sample(1:length(x),replace=TRUE)
      x.boot=x[rand.idx]
      y.boot=y[rand.idx]
      if(method == "spearman"){
        boot.sim[i] = cor(x.boot, y.boot, method="spearman", use="complete.obs")
      }else if(method == "pearson"){
        boot.sim[i] = cor(x.boot, y.boot, method="pearson",use="complete.obs")
      }else if(method == "bray"){
        boot.sim[i] = vegdist(rbind(x.boot,y.boot),method="bray")
      }else if(method == "kld"){
        boot.sim[i] = computeKld(x.boot,y.boot)
      }
    }
    boot.sim = na.omit(boot.sim)
    if(plot == T){
      hist(boot.sim,col=col2,add=T)
      abline(v=mean(boot.sim),col="red")
      legend(x="topleft", c("Permut","Boot"), bg="white",col=c(col1,col2),lty=rep(1,2),merge=T)
    }
    # if we got enough non-NA permutation and bootstrap values, compute p-value
    if(length(rand.sim) > round(N.rand/3) && length(boot.sim) > round(N.rand/3)){
      pval = pnorm(mean(rand.sim),mean=mean(boot.sim),sd=sd(boot.sim), lower.tail=lower.tail)
    }else{
      pval = 0.5
    }
  }else{
    # if we got enough non-NA permutation values, compute p-value
    if(length(rand.sim) > round(N.rand/3)){
      if (lower.tail) {
        pval = (sum(this.sim > rand.sim) / length(rand.sim))
      } else {
        pval = (sum(this.sim < rand.sim) / length(rand.sim))
      }
    }else{
      pval = 0.5
    }
  }
  # set missing value (from constant vector) to intermediate p-value (worst possible p-value in this context)
  if(is.na(pval)){
    pval = 0.5
  }
  # p-values are one-sided, so high p-values signal mutual exclusion and are converted into low ones
  if(pval > 0.5){
    pval = 1 - pval
  }
  if(verbose == T){
    print(paste("p-value =",pval))
    print(paste("original score",this.sim))
    print(paste("mean of null distrib",mean(rand.sim)))
    print(paste("sd of null distrib",sd(rand.sim)))
    if(permutandboot == T){
      print(paste("mean of boot distrib",mean(boot.sim)))
      print(paste("sd of boot distrib",sd(boot.sim)))
    }
  }
  pval
}
