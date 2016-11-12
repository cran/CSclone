#' @name Hamming
#' @title Hamming distance
#' @aliases Hamming
#' @description Hamming is computing the co-occurrence matrix.
#' @param k.matrix is matrix with the result-retaining of MCMC.
#' @details Hamming use the result-retaining of MCMC to compute the co-occurrence matrix. And we use the co-occurrence matrix to cluster which the criteria is MPEAR.
#' @return Hamming.matrix is the co-occurrence matrix and the idea is Hamming distance.
#' @return pro.matrix is the standardization of co-occurrence matrix.
#' @return data is data.frame with matrix. The purpose is constructing the heat map.
#' @author Peter Wu (peter123wu0@gmail.com)
#' @export
Hamming=function(k.matrix){
  cc=ncol(k.matrix)
  rr=nrow(k.matrix)
  Hamming.matrix=matrix(0,nrow=cc,ncol=cc)
  for(i in 1:(cc-1)){
    for(j in (i+1):cc){
      Hamming.matrix[i,j]=sum(k.matrix[,i]!=k.matrix[,j])      
    }
  }
  Hamming.matrix=Hamming.matrix+t(Hamming.matrix)
  pro.matrix=1-(Hamming.matrix/rr)
  Data=data.frame(x=rep(1:cc,each=cc),y=rep(1:cc,cc),value=as.vector(pro.matrix))
  return(list(Hamming.matrix=Hamming.matrix,pro.matrix=pro.matrix,data=Data))
}