#' @import lpSolve
#' @import mcclust
#' @import moments
#' @import DNAcopy
#' @import stats

#' @name CSclone
#' @title Clustering model which contains the CNAs mutation and SNVs mutation.
#' @aliases CSclone
#' @description CSclone is the main function of the package and a clustering model which contains the CNAs mutation and SNVs mutation.
#' @param somatic.id is a vector of s observiations and mark which loci has somatic mutation.
#' @param DNAcopy.object is a list and output of the DNAcopy package.
#' @param mcmc is a list giving the MCMC parameters. The list must include the following integers: nburn giving the number of burn-in scans, nskip giving the thinning interval, nsave giving the total number of scans to be saved, and ndisplay giving the number of saved scans to be displayed on screen (the function reports on the screen when every ndisplay iterations have been carried out) and default is list(nburn=5000,nsave=10000,nskip=1,ndisplay=1000).
#' @param y is a matrix giving the binomial data. The first column is B allele read depth and the second column is total read depth.
#' @param set is number that you can reproduce the simulation result and default is 100.
#' @param alpha is a number giving the concentration parameter of Dirichlet process and default is 1.
#' @param max.k is a number giving limit the maximum cluster number and default is NULL.
#' @param method is the agglomeration method to be used and the method is used to result of Hamming distance. The default is ward.D2.
#' @param prior is a vector has tow number which is the parameter of Beta distribution and the default is (1,1).
#' @details CSclone is a two steps model. The first step is clustering the non-CNAs somatic mutations by DPMM with Binomial distribution and the posterior distribution is Beta distribution. The second step classifies the CNAs mutation and the somatic mutation with CNAs mutation. There are an example and the dataset is small case in order to run quickly.
#' @return cluster is a vector giving the group result of the somatic mutation.
#' @return segment is a matrix. The first column is chromosome, the second column is starting position, the third column is end position, the fourth column is copy number, the fifth column is proportion of mutation, the sixth is number of loci, the seventh column is starting loci number, and the eighth column is end loci number.
#' @return mcmc is a list giving the MCMC parameters.
#' @return alpha is number giving the concentration parameter of Dirichlet process.
#' @return y is a matrix giving the binomial data.
#' @return group.prop is number giving the predict the proportion of group.
#' @return SNV is matrix giving the somatic mutation of every loci. The first column is B allele read depth, the second column is total read depth, the third column proportion of somatic mutation.
#' @author Peter Wu (peter123wu0@gmail.com)
#' @examples
#' mcmc=list(nburn=200,nsave=500,nskip=1,ndisplay=1000)
#' p=c(0.2,0.4,0.6,0.8)
#' chrs=rep(c(1:2),times=rep(500,2))
#' pos=sort(sample(size=1000,x=1:10^7))
#' pc=sample(rep(p,4*c(0,0,0.5,0.5)))
#' ps=sample(rep(p,50*c(0.3,0.3,0.2,0.2)))
#' x=simu.data(n.germline=1000,pc=pc,read=200,ps=ps,dis="Negative binomial",parameter=0.75)
#' snv.id=x$snv.id
#' y=x$y
#' row.names(y)=paste0(chrs,"_",pos)
#' logR=log(y[,2],base=2)-median(log(y[,2],base=2))
#' data=data.frame(chr=chrs,pos=pos,logR=logR)
#' rownames(data)=paste0("SNP",1:nrow(data))
#' CNA.object=CNA(data$logR,data$chr,data$pos,data.type="logratio",sampleid="test")
#' smoothed.CNA.object=smooth.CNA(CNA.object)
#' DNAcopy.object=segment(smoothed.CNA.object, undo.splits = "sdundo",undo.SD = 3, verbose = 1)
#' fit=CSclone(somatic.id=snv.id,DNAcopy.object=DNAcopy.object,mcmc=mcmc,y=y)
#' result=fit$cluster
#' tf_similar(real=ps,cluster=result)
#' 
#' @export
CSclone <- function(somatic.id,DNAcopy.object,mcmc=list(nburn=5000,nsave=10000,nskip=1,ndisplay=1000),y,set=100,alpha=1,max.k=NULL,method="ward.D2",prior=c(1,1)){
  STR <-median(y[,2])
  #set can fix random time
  set.seed(set)
  
  ###detect cna
  nc <-nrow(DNAcopy.object$output)
  x <-DNAcopy.object$output
  seg.ratio <-x$seg.mean
  seg.num <-x$num.mark
  seg.pos <-cumsum(x$num.mark)
  seg.id <-matrix(nrow=nc,ncol=2)
  for(i in 1:nc){
    if(i==1){
      seg.id[i,1] <-1
      seg.id[i,2] <-seg.pos[i]
    }else{
      seg.id[i,1] <-seg.pos[i-1]+1
      seg.id[i,2] <-seg.pos[i]      
    }  
  }
  seg.read <-sapply(1:nc,function(i){
    median(y[seg.id[i,1]:seg.id[i,2],2])
  })
  
  
  cns <-sapply(1:nc,function(i){
    if(seg.ratio[i]< (-0.1)){1
    }else if(seg.ratio[i]>=-0.1 & seg.ratio[i]<=0.1){2
    }else{3}
  })
  
  cna.id <-c()
  cna.id <-unlist(sapply(1:nc,function(i){
    c(cna.id,rep(cns[i],seg.id[i,2]-seg.id[i,1]+1))
  }))
  deletion.id <-which(cna.id==1)
  no.cna.id <-which(cna.id==2)
  amplification.id <-which(cna.id==3)
  cluster.id <-intersect(no.cna.id,somatic.id)
  s <-y[cluster.id,]
  y <-cbind(y,cna.id)
  
  ns <-nrow(s)
  cnt <-rep(2,ns)
  cnb1 <-rep(0,ns)
  cnb2 <-sapply(1:ns,function(i){
    if(s[i,1]/s[i,2]>0.5){2}else{1}
  })
  ###take the no cna and snv locus to cluster for DP
  ps <-sapply(1:ns,function(i){l.s(cnt=cnt[i],cnb1=cnb1[i],cnb2=cnb2[i],baf=s[i,1]/s[i,2])})
  
  ###mcmc
  run <-mcmc$nsave
  mc.id <-seq(from=(mcmc$nburn+1),to=run,by=mcmc$nskip)     
  #run is how many iteration
  #mc.id is which iteration retained
  
  new.p <-alpha*choose(s[,2],s[,1])*beta(s[,1]+prior[1],(s[,2]-s[,1])+prior[2])/beta(prior[1],prior[2])
  k.matrix <-matrix(ncol=ns,nrow=run)
  #k.matrix is the result for DP(MCMC)
  #nk is how many item in k group
  
  ###start DP
  #run=1
  k <-1
  kk <-1
  nk <-1
  for(i in 2:ns){
    pk <-c()
    for(m in 1:k){
      k.id <-which(kk==m)
      p.k <-median(ps[k.id])
      baf <-l.s(cnt=cnt[i],cnb1=cnb1[i],cnb2=cnb2[i],ps=p.k)
      pk[m] <-dbinom(s[i,1],s[i,2],baf)*nk[m] 
    }
    if(sum(pk)==0){
      pk[1] <-10^-100
    }
    prob <-c(pk,new.p[i])
    test <-which(rmultinom(1,size=1,prob=prob)==1)
    if(test==(k+1)){
      k <-k+1
      nk[k] <-1
      kk[i] <-k
    }else{
      nk[test] <-nk[test]+1
      kk[i] <-test
    }
  }
  k.matrix[1,] <-kk
  
  ###run>1
  for(r in 2:run){
    #Gibbs sampling for i locus then we need remove i locus imformation for group
    for(i in 1:ns){        
      pk <-c()
      nk[kk[i]] <-nk[kk[i]]-1
      kk[i] <-0     
      for(m in 1:k){           
        if(nk[m]==0){            
          pk[m] <-0
        }else{ 
          k.id <-which(kk==m)
          p.k <-median(ps[k.id])
          baf <-l.s(cnt=cnt[i],cnb1=cnb1[i],cnb2=cnb2[i],ps=p.k)
          pk[m] <-dbinom(s[i,1],s[i,2],baf)*nk[m]         
        }
      }
      if(sum(pk)==0){
        pk[1] <-10^-100
      }
      prob <-c(pk,new.p[i])
      test <-which(rmultinom(1,size=1,prob=pk)==1)        
      nk[test] <-nk[test]+1
      kk[i] <-test
    }
    
    k.matrix[r,] <-kk
    
    ###avoid having empty position for number of group
    if(length(unique(kk))!=max(kk)){
      nk <-c()
      for(j in 1:length(unique(kk))){
        kk[which(kk==sort(unique(kk))[j])] <-j
        nk[j] <-length(which(kk==j))
      }
    }
    k <-length(nk)   
    
    if(((r/mcmc$ndisplay)%% 1)==0){cat("MCMC ",r," of ",run,"\n")}
  }
  
  k.matrix <-k.matrix[mc.id,]
  fit.mm <-Hamming(k.matrix)
  mx <-adjust.mpear(fit.mm$pro.matrix,max.k=max.k)
  
  ###classify
  group.number <-length(unique(mx$cl))
  group.prop <-sapply(1:group.number,function(i){
    median(ps[mx$cl==i])
  })
  group.prop0 <-c(group.prop,group.prop,0)
  
  
  ###classify the proportion of CNAs mutation
  seg.pc <-c()
  seg.cna <-c()
  deletion.read <-c(group.prop*1+(1-group.prop)*2,group.prop*0+(1-group.prop)*2,2)*STR/2
  deletion.cna <-c(rep(1,group.number),rep(0,group.number),2)
  
  amplification.read <-c(group.prop*3+(1-group.prop)*2,group.prop*4+(1-group.prop)*2,2)*STR/2
  amplification.cna <-c(rep(3,group.number),rep(4,group.number),2)
  
  for(i in 1:length(cns)){
    if(cns[i]==1){
      seg.pc[i] <-group.prop0[which.min(abs(seg.read[i]-deletion.read))]
      seg.cna[i] <-deletion.cna[which.min(abs(seg.read[i]-deletion.read))]
    }else if(cns[i]==2){
      seg.pc[i] <-0
      seg.cna[i] <--1
    }else if(cns[i]==3){
      seg.pc[i] <-group.prop0[which.min(abs(seg.read[i]-amplification.read))]
      seg.cna[i] <-amplification.cna[which.min(abs(seg.read[i]-amplification.read))]
    }
  }
  
  cna <-unlist(sapply(1:nc,function(i){
    rep(seg.cna[i],seg.num[i])
  }))
  
  cna.pc <-unlist(sapply(1:nc,function(i){
    rep(seg.pc[i],seg.num[i])
  }))
  
  ###classify the proportion of SNVs mutation which has CNAs mutation  
  #cnt=0
  id.0 <-intersect(which(cna==0),somatic.id)
  if(length(id.0)==1){
    read.0 <-matrix(nrow=1,ncol=2,data=y[id.0,-3])
    colnames(read.0) <-c("B.read","total.read")
  }else{read.0 <-y[id.0,-3]}
  
  cluster.0 <-c()
  if(length(id.0)>0){
    for(i in 1:length(id.0)){
      like <-c()
      for(j in 1:group.number){
        #AA->AB->0
        if(group.prop[j]>cna.pc[id.0[i]]){
          baf1 <-l.s(cnt=0,cnb1=0,cnb2=1,pc=cna.pc[id.0[i]],ps=group.prop[j])
        }else{baf1 <-0}
        #AA->BB->0
        if(group.prop[j]>cna.pc[id.0[i]]){
          baf2 <-l.s(cnt=0,cnb1=0,cnb2=2,pc=cna.pc[id.0[i]],ps=group.prop[j])
        }else{baf2 <-0}
        baf <-c(baf1,baf2)
        like <-c(like,dbinom(read.0[i,1],read.0[i,2],baf))
      }
      cluster.0[i] <-rep(1:group.number,each=2)[which.max(like)]
    }
  }
  
  #cnt=1
  id.1 <-intersect(which(cna==1),somatic.id)
  if(length(id.1)==1){
    read.1 <-matrix(nrow=1,ncol=2,data=y[id.1,-3])
    colnames(read.1) <-c("B.read","total.read")
  }else{read.1 <-y[id.1,-3]}
  
  cluster.1 <-c()
  if(length(id.1)>0){
    for(i in 1:length(id.1)){
      like <-c()
      for(j in 1:group.number){
        #AA->AB->A
        if(group.prop[j]>cna.pc[id.1[i]]){
          baf1 <-l.s(cnt=1,cnb1=0,cnb2=1,pc=cna.pc[id.1[i]],ps=group.prop[j])
        }else{baf1 <-0}
        #AA->AB->B
        if(group.prop[j]>cna.pc[id.1[i]]){
          baf2 <-l.s(cnt=1,cnb1=1,cnb2=1,pc=cna.pc[id.1[i]],ps=group.prop[j])
        }else{baf2 <-0}
        #AA->BB->B
        if(group.prop[j]>cna.pc[id.1[i]]){
          baf3 <-l.s(cnt=1,cnb1=1,cnb2=2,pc=cna.pc[id.1[i]],ps=group.prop[j])
        }else{baf3 <-0}
        #AA->A->B
        if(group.prop[j]<=cna.pc[id.1[i]]){
          baf4 <-l.s(cnt=1,cnb1=0,cnb2=1,pc=cna.pc[id.1[i]],ps=group.prop[j])
        }else{baf4 <-0}    
        baf <-c(baf1,baf2,baf3,baf4)
        like <-c(like,dbinom(read.1[i,1],read.1[i,2],baf))
      }
      cluster.1[i] <-rep(1:group.number,each=4)[which.max(like)]
    }
  }
  
  #cnt=2
  id.2 <-intersect(which(cna==2),somatic.id)
  if(length(id.2)==1){
    read.2 <-matrix(nrow=1,ncol=2,data=y[id.2,-3])
    colnames(read.2) <-c("B.read","total.read")
  }else{read.2 <-y[id.2,-3]}
  
  cluster.2 <-c()
  if(length(id.2)>0){
    for(i in 1:length(id.2)){
      like <-c()
      for(j in 1:group.number){
        #AA->AB
        baf1 <-l.s(cnt=2,cnb1=1,cnb2=0,ps=group.prop[j])
        #AA->BB
        baf2 <-l.s(cnt=2,cnb1=2,cnb2=0,ps=group.prop[j])
        baf <-c(baf1,baf2)
        like <-c(like,dbinom(read.2[i,1],read.2[i,2],baf))
      }
      cluster.2[i] <-rep(1:group.number,each=2)[which.max(like)]
    }
  }
  
  #cnt=3
  id.3 <-intersect(which(cna==3),somatic.id)
  if(length(id.3)==1){
    read.3 <-matrix(nrow=1,ncol=2,data=y[id.3,-3])
    colnames(read.3) <-c("B.read","total.read")
  }else{read.3 <-y[id.3,-3]}
  
  cluster.3 <-c()
  if(length(id.3)>0){
    for(i in 1:length(id.3)){
      like <-c()
      for(j in 1:group.number){
        #AA->AB->AAB
        if(group.prop[j]>cna.pc[id.3[i]]){
          baf1 <-l.s(cnt=3,cnb1=1,cnb2=1,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf1 <-0}
        #AA->AB->ABB
        if(group.prop[j]>cna.pc[id.3[i]]){
          baf2 <-l.s(cnt=3,cnb1=2,cnb2=1,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf2 <-0}
        #AA->BB->BBB
        if(group.prop[j]>cna.pc[id.3[i]]){
          baf3 <-l.s(cnt=3,cnb1=3,cnb2=2,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf3 <-0}
        #AA->AAA->AAB
        if(group.prop[j]<=cna.pc[id.3[i]]){
          baf4 <-l.s(cnt=3,cnb1=0,cnb2=1,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf4 <-0}
        #AA->AAA->ABB
        if(group.prop[j]<=cna.pc[id.3[i]]){
          baf5 <-l.s(cnt=3,cnb1=0,cnb2=2,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf5 <-0}
        #AA->AAA->BBB
        if(group.prop[j]<=cna.pc[id.3[i]]){
          baf6 <-l.s(cnt=3,cnb1=0,cnb2=3,pc=cna.pc[id.3[i]],ps=group.prop[j])
        }else{baf6 <-0}
        baf <-c(baf1,baf2,baf3,baf4,baf5,baf6)
        like <-c(like,dbinom(read.3[i,1],read.3[i,2],baf))
      }
      cluster.3[i] <-rep(1:group.number,each=6)[which.max(like)]
    }
  }
  
  #cnt=4
  id.4 <-intersect(which(cna==4),somatic.id)
  if(length(id.4)==1){
    read.4 <-matrix(nrow=1,ncol=2,data=y[id.4,-3])
    colnames(read.4) <-c("B.read","total.read")
  }else{read.4 <-y[id.4,-3]}
  
  cluster.4 <-c()
  if(length(id.4)>0){
    for(i in 1:length(id.4)){
      like <-c()
      for(j in 1:group.number){
        #AA->AB->AAAB
        if(group.prop[j]>cna.pc[id.4[i]]){
          baf1 <-l.s(cnt=4,cnb1=1,cnb2=1,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf1 <-0}
        #AA->AB->ABBB
        if(group.prop[j]>cna.pc[id.4[i]]){
          baf2 <-l.s(cnt=4,cnb1=3,cnb2=1,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf2 <-0}
        #AA->BB->BBBB
        if(group.prop[j]>cna.pc[id.4[i]]){
          baf3 <-l.s(cnt=4,cnb1=4,cnb2=2,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf3 <-0}
        #AA->AAAA->AAAB
        if(group.prop[j]<=cna.pc[id.4[i]]){
          baf4 <-l.s(cnt=4,cnb1=0,cnb2=1,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf4 <-0}
        #AA->AAAA->AABB
        if(group.prop[j]<=cna.pc[id.4[i]]){
          baf5 <-l.s(cnt=4,cnb1=0,cnb2=2,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf5 <-0}
        #AA->AAAA->ABBB
        if(group.prop[j]<=cna.pc[id.4[i]]){
          baf6 <-l.s(cnt=4,cnb1=0,cnb2=3,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf6 <-0}
        #AA->AAAA->BBBB
        if(group.prop[j]<=cna.pc[id.4[i]]){
          baf7 <-l.s(cnt=4,cnb1=0,cnb2=4,pc=cna.pc[id.4[i]],ps=group.prop[j])
        }else{baf7 <-0}
        baf <-c(baf1,baf2,baf3,baf4,baf5,baf6,baf7)
        like <-c(like,dbinom(read.4[i,1],read.4[i,2],baf))
      }
      cluster.4[i] <-rep(1:group.number,each=7)[which.max(like)]
    }
  }  
  
  #summarize the clustering result of SNVs mutation
  final.cluster <-c()
  for(i in 1:ns){final.cluster[which(somatic.id==cluster.id[i])] <-mx$cl[i]}  
  for(i in 1:length(id.0)){final.cluster[which(somatic.id==id.0[i])] <-cluster.0[i]}
  for(i in 1:length(id.1)){final.cluster[which(somatic.id==id.1[i])] <-cluster.1[i]}
  for(i in 1:length(id.2)){final.cluster[which(somatic.id==id.2[i])] <-cluster.2[i]}  
  for(i in 1:length(id.3)){final.cluster[which(somatic.id==id.3[i])] <-cluster.3[i]}
  for(i in 1:length(id.4)){final.cluster[which(somatic.id==id.4[i])] <-cluster.4[i]} 
  seg <-data.frame(chr=x[,2],loc.start=x[,3],loc.end=x[,4],cna=seg.cna,proportion=seg.pc,number=seg.num,from=seg.id[,1],to=seg.id[,2])
  ss <-seg[1,]
  a <-1
  i <-1
  while(i<nrow(seg)){
    b <-1
    if(seg[i,1]!=seg[i+1,1] | seg[i,4]!=seg[i+1,4] | seg[i,5]!=seg[i+1,5]){
      ss[a,] <-seg[i,]
      a <-a+1
    }else{
      while(seg[i,1]==seg[i+b,1] & seg[i,4]==seg[i+b,4] & seg[i,5]==seg[i+b,5] & (i+b)<=nrow(seg)){
        ss[a,] <-seg[i,]
        ss[a,3] <-seg[i+b,3]
        ss[a,6] <-sum(seg[i:i+b,6])
        ss[a,8] <-seg[i+b,8]
        b <-b+1
      }
      a <-a+1
    }
    i <-i+b
  }
  if(b==1){ss[a,] <-seg[nrow(seg),]}
  prop=sapply(1:length(somatic.id),function(i){group.prop[final.cluster[i]]})
  SNV <-cbind(y[somatic.id,-3],prop)
  return(list(cluster=final.cluster,segment=ss,mcmc=mcmc,alpha=alpha,y=y,group.prop=group.prop,SNV=SNV))
}