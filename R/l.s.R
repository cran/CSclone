#' @name l.s
#' @title Linear relation for SNVs
#' @aliases l.s
#' @description l.s is the linear relation between B allele frequence and the proportion of SNVs mutation(PS).
#' @param cnt is number with the total copy number after the CNAs mutation.
#' @param cnb1 is number with the B copy number of two types of periods. The first period, the proportion of CNAs mutation is more than the proportion of SNVs mutation then cnb1 is B copy number after SNVs mutation. The second period, the proportion of CNAs mutation is less than the proportion of SNVs mutation then cnb1 is B copy number between the SNVs mutation and CNAs mutation.
#' @param cnb2 is number with the B copy number of two types of periods. The first period, the proportion of CNAs mutation is less than the proportion of SNVs mutation then cnb1 is B copy number after CNAs mutation. The second period, the proportion of CNAs mutation is more than the proportion of SNVs mutation then cnb1 is B copy number between the CNAs mutation and SNVs mutation.
#' @param pc is number with the proportion of CNAs mutation and the default is 0.
#' @param ps is number with the proportion of SNVs mutation and the default is NULL.
#' @param baf is number with the B allele frequence and the default is NULL.
#' @details l.s is bidirection function. The first function is given the proportion of SNVs mutation(ps) to predict the B allele frequence(baf). The second function is given the B allele frequence(baf) to predict the proportion of SNVs mutation(ps).
#' @return baf is number with B allele frequence(baf) if the input is given the proportion of SNVs mutation(ps).
#' @return ps is number with the proportion of SNVs mutation(ps) if the input is given the B allele frequence(baf).
#' @author Peter Wu (peter123wu0@gmail.com)
#' @export
l.s=function(cnt,cnb1,cnb2,pc=0,ps=NULL,baf=NULL){
  # linear relation between baf and ps condition of pc
  if(is.null(baf)){
    if(pc==0){
      baf=ps*cnb2/2
    }else if(pc>=ps){
      baf=ps*cnb2/(2*(1-pc)+pc*cnt)
    }else{
      baf=((ps-pc)*cnb2+pc*cnb1)/(2*(1-pc)+pc*cnt)
    }
    if(baf>1)baf=1
    return(baf)
  }
  if(is.null(ps)){
    if(pc==0){
      ps=baf*2/cnb2
    }else{
      ps=baf*((2*(1-pc)+pc*(cnt))/cnb2)
      if(ps>pc){
        ps=pc+(baf*(2*(1-pc)+pc*cnt)-pc*cnb1)/cnb2
      }
    }
    return(ps)
  }
}