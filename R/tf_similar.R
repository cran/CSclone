#' @name tf_similar
#' @title True or false of similar
#' @aliases tf_similar
#' @description tf_similar is comparing the similarity of different group.
#' @param real is vector with the real group number. If we don't have the real group number, we also can input the result of some clustering method.
#' @param cluster is vector with the result of clustering.
#' @details tf_similar is comparing the similarity of different group. So the main purpose is checking the quality of clustering method. We show the rand index(RI), weight rand index(WRI), and adjust rand index(ARI) to compare. In order to compute the three kinds of index, we count the TT, TF, FT, FF, and the relation weight. The details of the function of the function  are in paper of Tumor subclones detection with Dirichlet Process Mixture Model.
#' @return result is data.frame with 10 columns which contains the TT, TF, FT, FF, rand index(RI), weight rand index(WRI), adjust rand index(ARI), the weighted TT, the weighted TF, and the weight FT.
#' @return weight is data.frame with 3 columns which contains the weight about the TT, TF, and FT.
#' @author Peter Wu (peter123wu0@gmail.com)
#' @export
tf_similar=function(real,cluster){ 
  if(length(real)!=length(cluster)){
    stop("the length need to be equal")
  }
  N=length(real)
  TT=sum(choose(table(real,cluster),2))
  TF=sum(choose(table(real),2))-TT
  FT=sum(choose(table(cluster),2))-TT
  FF=choose(N,2)-TF-FT-TT
  
  e_TT=sum(choose(table(real),2))*sum(choose(table(cluster),2))/choose(N,2)
  e_TF=sum(choose(table(real),2))-e_TT
  e_FT=sum(choose(table(cluster),2))-e_TT
  ee=(1/e_TT)+(1/e_TF)+(1/e_FT)
  w_TT=(1/e_TT)/ee
  w_TF=(1/e_TF)/ee
  w_FT=(1/e_FT)/ee
  
  #RI is Rand index
  RI=(TT+FF)/choose(N,2)  
  
  #WRI is weighted Rand index
  WRI=w_TT*TT/(w_TT*TT+w_TF*TF+w_FT*FT)
  
  #ARI is adjust Rand index
  ARI=(TT-e_TT)/(((TT+TF)+(TT+FT))/2-e_TT)
  
  result=as.data.frame(matrix(c(RI,TT,TF,FT,FF,ARI,WRI,w_TT*TT,w_TF*TF,w_FT*FT),nrow=length(TT),ncol=10))
  weight=as.data.frame(matrix(c(w_TT,w_TF,w_FT),nrow=length(TT),ncol=3))
  colnames(result)=c("RI","TT","TF","FT","FF","ARI","WRI","w.TT*TT","w.TF*TF","w.FT*FT")
  colnames(weight)=c("w_TT","w_TF","w_FT")
  return(list(result=result,weight=weight))
}