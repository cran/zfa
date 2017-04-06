zooming.fast<-function(data,y,obj,test=c("wtest","SKAT","SKATO","burden")) {
  lower<-1
  upper<-ncol(data)
  mid<-ncol(data)
  if (test=="wtest") {
    pvalue<-wtest(data[,c(lower:upper)],y)
  } else if (test=="SKAT") {
    pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower:upper)]),obj)$p.value
  } else if (test=="SKATO") {
    pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower:upper)]),obj,method="SKATO")$p.value
  } else if (test=="burden") {
    pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower:upper)]),obj,method="Burden")$p.value
  }
  pvalue.all<-pvalue
  upper.all<-upper
  lower.all<-lower
  while (mid>2){
    mid<-mid/2
    if (test=="wtest") {
      pvalue.left<-wtest(data[,c(lower:(lower+mid-1))],y)
      pvalue.right<-wtest(data[,c((upper-mid+1):upper)],y)
    } else if (test=="SKAT") {
      pvalue.left<-(SKAT::SKATBinary(as.matrix(data[,c(lower:(lower+mid-1))]),obj)$p.value)
      pvalue.right<-(SKAT::SKATBinary(as.matrix(data[,c((upper-mid+1):upper)]),obj)$p.value)
    } else if (test=="SKATO") {
      pvalue.left<-(SKAT::SKATBinary(as.matrix(data[,c(lower:(lower+mid-1))]),obj,method="SKATO")$p.value)
      pvalue.right<-(SKAT::SKATBinary(as.matrix(data[,c((upper-mid+1):upper)]),obj,method="SKATO")$p.value)
    } else if (test=="burden") {
      pvalue.left<-(SKAT::SKATBinary(as.matrix(data[,c(lower:(lower+mid-1))]),obj,method="Burden")$p.value)
      pvalue.right<-(SKAT::SKATBinary(as.matrix(data[,c((upper-mid+1):upper)]),obj,method="Burden")$p.value)
    }
    pvalue<-min(pvalue.left,pvalue.right)
    if (pvalue==pvalue.left){
      upper<-upper-mid
    } else if (pvalue==pvalue.right){
      lower<-lower+mid
    }
    pvalue<-ifelse(pvalue>1,1,pvalue)
    pvalue.all<-c(pvalue.all,pvalue)
    upper.all<-c(upper.all,upper)
    lower.all<-c(lower.all,lower)
  }
  pvalue.final<-min(pvalue.all)
  upper.final<-upper.all[which(pvalue.all==pvalue.final)[1]]
  lower.final<-lower.all[which(pvalue.all==pvalue.final)[1]]
  pvalue.final<-pvalue.final
  pvalue.final<-ifelse(pvalue.final>1,1,pvalue.final)
  result<-cbind(lower.final,upper.final,pvalue.final)
  colnames(result)<-c("lower","upper","corrected.pvalue")
  rm(pvalue,pvalue.all,pvalue.final,upper,upper.all,upper.final,lower,lower.all,lower.final,data,y)
  return(result)
}
