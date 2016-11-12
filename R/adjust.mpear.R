#' @name adjust.mpear
#' @title Adjust MPEAR
#' @aliases adjust.mpear
#' @description adjust.mpear is modified the mcclust package to suit our model.
#' @param psm is posterior similarity matrix with entries between 0 and 1 and 1's on the diagonals which is distance matrix.
#' @param method is clustering method. We offer many kinds of method, like as kmeans and hierarchical clustering. The default is ward.D2 of hierarchical clustering.
#' @param max.k is number with the user has the professional judgment to consider the max group and the default is NULL.
#' @details adjust.mpear is modified the mcclust package to suit our model. And the function use distance matrix by criteria of MPEAR to cluster the data.
#' @return cl is vector with the clustering result.
#' @return value is number with MPEAR.
#' @return method is clustering method. 
#' @author Peter Wu (peter123wu0@gmail.com)
#' @references Fritsch, A. and Ickstadt, K. Improved criteria for clustering based on the posterior similarity matrix. Bayesian analysis 2009;4(2):367-391.
#' @export
adjust.mpear=function(psm, method = "ward.D2" , max.k = NULL){
  if(any(psm != t(psm)) | any(psm > 1) | any(psm < 0) | sum(diag(psm)) != nrow(psm)){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")
  }
  
  if(is.null(max.k))max.k=10
  k=1
  a=0
  while(a<0.99 & k<=max.k){
    a=kmeans(1-psm,k)$betweenss/kmeans(1-psm,k)$totss
    k=k+1
  }
  k=k-1
  
  if(sum(method == c("kmeans","single","average","complete","ward.D2")) == 0){
    stop("pairwise cluster method is not right option, consider kmeans or HC")
  }
  
  if(method == "kmeans"){
    cls.kmeans = t(apply(matrix(1:k), 1, function(x){ kmeans(1-psm,x)$cluster })) 
    pears.kmeans = pear(cls.kmeans, psm)
    val.kmeans = max(pears.kmeans)
    cl.kmeans = cls.kmeans[which.max(pears.kmeans), ]
    return(list(cl = cl.kmeans, value = val.kmeans, method = "kmeans"))
  }else{
    hc = hclust(as.dist(1 - psm), method = method)
    cls.hc = t(apply(matrix(1:k), 1, function(x){ cutree(hc, k = x) })) 
    pears.hc = pear(cls.hc, psm)
    val.hc = max(pears.hc)
    cl.hc = cls.hc[which.max(pears.hc), ]
    return(list(cl = cl.hc, value = val.hc, method = method))
  }
}