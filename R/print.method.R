#' @export

print.zfa<-function(x, ...){
  cat("No.regions:\n")
  print(x$n.region)
  cat("\nNo.rare:\n")
  print(x$n.rare)
  cat("\nNo.common:\n")
  print(x$n.common)
  cat("\nResults:\n")
  print(x$results)
  cat("Bon.sig.level:\n")
  print(x$Bon.sig.level)
}

