#' @name simu.data
#' @title Simulate data
#' @aliases simu.data
#' @description simu.data is generating the data which simulates the locus. 
#' @param n.germline is number with the total number of generating data which contains non-CNAs and non-SNVs. The default is 10000.
#' @param segment.length is number with the number of loci of every segment and the number is common divisor of n.germline. The default is 100.
#' @param seed is number that you can reproduce the simulation result and default is NULL.
#' @param pc is vector with the proportion of CNAs mutation of every segment.
#' @param read is number with the standard total read depth.
#' @param ps is vector with the proportion of SNVs mutation.
#' @param dis is distribution with generating the data and the default is Normal distribution. We provide Poisson and Negative binomial distribution to choose.
#' @param parameter is number and default is NULL. If dis is Normal and parameter is not NULL, then the standard deviation is parameter x mean. If dis is Normal and parameter is NULL, then the standard deviation is mean. If dis is Negative binomial and parameter is NULL, then the parameter is second parameter. If dis is Negative binomial and parameter NULL, then the second parameter is 0.5.
#' @details The simu.data can consider not only SNVs mutation but also CNAs mutation and generate the data which contains germline data and mutation data.
#' @return y is n x 2 matrix. The first column is the B allele read depth and the second column is the total read depth.
#' @return snv.id is vector which denotes which loci with SNVs mutation.
#' @return cnt.id is A x B matrix which denotes which loci with CNAs mutation.. The A is segment.length and the B is the number of segment with CNAs mutation.
#' @return cnv is n x 4 matrix. The first column is the total copy number after CNAs mutation, the second column is the cnb1, the third column is the cnb2, and the fourth column is the proportion of CNAs mutation. The cnb1 and the cnb2 has detailed explanation at l.s function.
#' @author Peter Wu (peter123wu0@gmail.com)
#' @export
simu.data=function(n.germline=10000,segment.length=100,seed=NULL,pc,read=100,ps,dis="Normal",parameter=NULL){
  n.seg=length(pc)
  n.locus=length(ps)
  if(is.null(seed)!=T){set.seed(seed)}else{set.seed(100)}
  cnt.number=sapply(1:n.seg,function(i){
    a=which(rmultinom(n.seg,1,c(0.5,0.4,0.1))[,i]==1)
    if(a>=2){a+1}else{a}
  })
  
  cnt.group=sort(sample(1:floor(n.germline/segment.length),n.seg))
  
  cnt.id=sapply(1:n.seg,function(i){
    ((cnt.group[i]-1)*segment.length+1):(cnt.group[i]*segment.length)
  })
  
  cnt=rep(2,n.germline)
  cnt.pc=rep(0,n.germline)
  for(i in 1:n.seg){
    cnt[cnt.id[,i]]=cnt.number[i]
    cnt.pc[cnt.id[,i]]=pc[i]
  }
  
  x=read*((2*(1-cnt.pc)+cnt*cnt.pc)/2)
  
  if(sum(dis == c("Normal","Poisson","Negative binomial")) == 0){
    stop("the distribution only provide Normal, Poisson, Negative binomial")
  }
  if(dis=="Normal"){
    total.read=sapply(1:n.germline,function(i){
      if(is.null(parameter)){
        round(rnorm(1,x[i],sqrt(x[i])))
      }else{
        round(rnorm(1,x[i],parameter*sqrt(x[i])))
      }
    })
  }
  
  if(dis=="Poisson"){
    total.read=sapply(1:n.germline,function(i){
      rpois(1,x[i])
    })
  }
  
  if(dis=="Negative binomial"){
    total.read=sapply(1:n.germline,function(i){
      if(is.null(parameter)){
        rnbinom(1,x[i],0.5)
      }else{
        rnbinom(1,parameter*x[i]/(1-parameter),parameter)
      }
    })
  }
  
  snv.id=sort(sample(n.germline,n.locus))
  locus.ps=rep(0,n.germline)
  locus.ps[snv.id]=ps
  
  cnb1=cnb2=B.read=rep(0,n.germline)
  for(i in snv.id){
    if(cnt.pc[i]>=locus.ps[i]){
      if(cnt[i]==0){
      }else{
        if(locus.ps[i]<=0.5){
          cnb2[i]=1
        }else{
          cnb2[i]=sapply(cnt[i],function(x)sample(1:x,1))
        }
      }
    }else{
      if(locus.ps[i]<=0.5){
        cnb2[i]=1
      }else{
        cnb2[i]=sample(size=1,x=1:2)
      }
      if(cnt[i]<2){
        cnb1[i]=cnt[i]
      }else{
        if(cnb2[i]==1){
          cnb1[i]=1
        }else{
          cnb1[i]=cnt[i] 
        }
      }
    }
    baf=l.s(cnt=cnt[i],cnb1=cnb1[i],cnb2=cnb2[i],ps=locus.ps[i],pc=cnt.pc[i])
    B.read[i]=rbinom(n=1,size=total.read[i],prob=baf)
  }
  
  cnv=cbind(cnt,cnb1,cnb2,cnt.pc)
  y=cbind(B.read,total.read)  
  return(list(y=y,snv.id=snv.id,cnt.id=cnt.id,cnv=cnv))
}