#'@useDynLib zfa, .registration = TRUE

table.e1<-function(x,y){
  .Call("table_e1",x,y)
}

