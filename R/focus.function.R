focusing<-function(result,data,y,fast.path,filter.pval,bin,i,obj,test=c("wtest","SKAT","SKATO","burden")) {
  filter.pval<-ifelse(is.null(filter.pval),1,filter.pval)
  if (result[,"corrected.pvalue"]<=filter.pval){
    d<-result[,"upper"]-result[,"lower"]+1
    if (d==2) {
      inward<-0
      outward<-2*d
    } else {
      inward<-d/4
      outward<-d/2
    }
    pvalue.left<-NULL
    lb.all<-NULL
    for (j in (-outward):(inward)) {
      lb.new<-result[,"lower"]+j
      if (lb.new>0) {
        lb.all<-rbind(lb.all,lb.new)
        if (test=="wtest") {
          pvalue.left<-rbind(pvalue.left,wtest(data[,c(lb.new:result[,"upper"])],y))
        } else if (test=="SKAT") {
          pvalue.left<-rbind(pvalue.left,SKAT::SKATBinary(as.matrix(data[,c(lb.new:result[,"upper"])]),obj)$p.value)
        } else if (test=="SKATO") {
          pvalue.left<-rbind(pvalue.left,SKAT::SKATBinary(as.matrix(data[,c(lb.new:result[,"upper"])]),obj,method="SKATO")$p.value)
        } else if (test=="burden") {
          pvalue.left<-rbind(pvalue.left,SKAT::SKATBinary(as.matrix(data[,c(lb.new:result[,"upper"])]),obj,method="Burden")$p.value)
        }
      }
    }
    lb.final<-lb.all[which(pvalue.left==min(pvalue.left))][1]
    rm(pvalue.left,lb.new)
    pvalue.right<-NULL
    ub.all<-NULL
    for (j in (-inward):(outward)) {
      ub.new<-result[,"upper"]+j
      if (ub.new<(ncol(data)+1)) {
        ub.all<-rbind(ub.all,ub.new)
        if (test=="wtest") {
          pvalue.right<-rbind(pvalue.right,wtest(data[,c(lb.final:ub.new)],y))
        } else if (test=="SKAT") {
          pvalue.right<-rbind(pvalue.right,SKAT::SKATBinary(as.matrix(data[,c(lb.final:ub.new)]),obj)$p.value)
        } else if (test=="SKATO") {
          pvalue.right<-rbind(pvalue.right,SKAT::SKATBinary(as.matrix(data[,c(lb.final:ub.new)]),obj,method="SKATO")$p.value)
        } else if (test=="burden") {
          pvalue.right<-rbind(pvalue.right,SKAT::SKATBinary(as.matrix(data[,c(lb.final:ub.new)]),obj,method="Burden")$p.value)        }
      }
    }
    ub.final<-ub.all[which(pvalue.right==min(pvalue.right))][1]
    k<-length(c(lb.all,ub.all))-2
    if (fast.path) {
      pvalue.right.corrected<-pvalue.right*(bin/2+log2(bin)-1+k)
    } else {
      pvalue.right.corrected<-pvalue.right*(bin-1+k)
    }
    pvalue.right.corrected<-ifelse(pvalue.right.corrected>1,1,pvalue.right.corrected)
    result<-cbind(lb.final+bin*(i-1),ub.final+bin*(i-1),ub.final-lb.final+1,min(pvalue.right.corrected)[1])
    colnames(result)<-c("lower.bound","upper.bound","opt.region.size","corrected.pvalue")
    rm(pvalue.right,pvalue.right.corrected,ub.all,ub.new)
  } else {
    if (fast.path) {
      result[,"corrected.pvalue"]<-result[,"corrected.pvalue"]*(bin/2+log2(bin)-1)
      result[,"corrected.pvalue"]<-ifelse(result[,"corrected.pvalue"]>1,1,result[,"corrected.pvalue"])
    } else {
      result[,"corrected.pvalue"]<-result[,"corrected.pvalue"]/(bin/(result[,"upper"]-result[,"lower"]+1))*(bin-1)
      result[,"corrected.pvalue"]<-ifelse(result[,"corrected.pvalue"]>1,1,result[,"corrected.pvalue"])
    }
    result<-t(c(result[,"lower"]+bin*(i-1),result[,"upper"]+bin*(i-1),result[,"upper"]-result[,"lower"]+1,result[,"corrected.pvalue"]))
    colnames(result)<-c("lower.bound","upper.bound","opt.region.size","corrected.pvalue")
  }
  return(result)
}
