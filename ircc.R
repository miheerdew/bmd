ircc.kmeans <- function (X, Y, k, ...) {
  
  
}

ircc <- function (X, Y, nbmds, ..., method = "kmeans") {
  
  if (method == "kmeans") {
    
    res <- ircc.kmeans(X, Y, k = nbmds, ...)
    
  }
  
  return(res)
  
}