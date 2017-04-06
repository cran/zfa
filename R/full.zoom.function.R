zooming.full<-function(data,y,obj,test=c("wtest","SKAT","SKATO","burden")) {
  P<-ncol(data)
  R<-floor(log2(P))
  r<-c(1:R)
  nr<-c(2^(r-1))
  d<-P/nr
  result.all<-NULL
  for (j in 1:R) {
    d1<-d[j]
    nr1<-nr[j]
    c<-c(1:nr1)
    lower.all<-d1*(c-1)+1
    upper.all<-d1*c
    pvalue.all<-NULL
    for (k in 1:nr1){
      if (test=="wtest") {
        pvalue<-wtest(data[,c(lower.all[k]:upper.all[k])],y)*nr1
      } else if (test=="SKAT") {
        pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower.all[k]:upper.all[k])]),obj)$p.value*nr1
      } else if (test=="SKATO") {
        pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower.all[k]:upper.all[k])]),obj,method="SKATO")$p.value*nr1
      } else if (test=="burden") {
        pvalue<-SKAT::SKATBinary(as.matrix(data[,c(lower.all[k]:upper.all[k])]),obj,method="Burden")$p.value*nr1
      }
      pvalue<-ifelse(pvalue>1,1,pvalue)
      pvalue.all<-c(pvalue.all,pvalue)
    }
    result<-cbind(lower.all,upper.all,pvalue.all)
    result.all<-rbind(result.all,result)
  }
  pvalue.final<-result.all[,3]
  pvalue.final<-ifelse(pvalue.final>1,1,pvalue.final)
  result.all<-cbind(result.all[,1:2],pvalue.final)
  colnames(result.all)<-c("lower","upper","corrected.pvalue")
  order.pvalue<-order(result.all[,"corrected.pvalue"],decreasing=F)
  min.pvalue<-t(result.all[order.pvalue[1],])
  rm(order.pvalue,result.all,result)
  return(min.pvalue)
}
