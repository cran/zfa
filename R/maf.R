maf<-function(data){
  p<-colMeans(data)/2
  p[p>0.5]<-1-p[p>0.5]
  result.maf<-p
  return(result.maf)
}
