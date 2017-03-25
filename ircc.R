# Examples:
#
# Specifiying the centers within each bimodule
# res <- ircc(sim$X, sim$Y, 4, centers.X=c(1,101,201,301), centers.Y=c(1,101,201,301))
#
# By specifying just the number of bimodules.
# res <- ircc(sim$X, sim$Y, 4)

ircc.kmeans <- function(X, Y, k, centers.X=NULL, centers.Y=NULL, ...) {
  cormat <- cor(X,Y)
  d <- ncol(X)

  as.partition <- function(partition_map){
    lapply(1:max(partition_map), function(i) which(partition_map == i))
  }


  if (is.null(centers.X)) {
    centers <- k
  } else {
    centers <- cormat[centers.X,]
  }
  
  CX <- as.partition(kmeans(cormat, centers, ...)$cluster)
  
  if (is.null(centers.Y)) {
    centers <- k
  } else {
    centers <- t(cormat)[centers.Y,]
  }
  CY <- as.partition(kmeans(t(cormat), centers, ...)$cluster)

  mean_cor <- function(i,j){
    mean(cormat[CX[[i]], CY[[j]]])
  }

  #Calculate mean-correlation between CX and CY.
  M <- outer(1:k, 1:k, function(x,y) mapply(mean_cor, x, y, SIMPLIFY = TRUE))

  #image(M)

  comms <- rep(list(NULL), k)
  X_ind <- Y_ind <- 1:k #Cluster indices
  count <- 1

  while (count <= k) {
    argmax <- arrayInd(which(M == max(M[X_ind, Y_ind])), .dim = c(k, k))[1,]
    X_ind <- setdiff(X_ind, argmax[1])
    Y_ind <- setdiff(Y_ind, argmax[2])
    comms[[count]] <- c(CX[[argmax[1]]], CY[[argmax[2]]] + d)
    count <- count + 1
  }
  return(comms)

}

ircc <- function(X, Y, nbmds, ..., method = "kmeans") {
  
  if (method == "kmeans") {
    
    res <- ircc.kmeans(X, Y, k = nbmds, ...)
    
  }
  
  return(res)
  
}