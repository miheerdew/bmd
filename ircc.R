# Examples:
#
# K-means
# ---------
# Specifiying the centers within each bimodule
# res <- ircc(sim$X, sim$Y, 4, centers.X=c(1,101,201,301), centers.Y=c(1,101,201,301))
#
# By specifying just the number of bimodules.
# res <- ircc(sim$X, sim$Y, 4)
#
#
# Hclust
# -------

clustering_approach_to_bmd <- function(X, Y, k, clustering_algo,
                                      Xargs=NULL,
                                      Yargs=NULL){
  cormat <- cor(X,Y)
  n <- nrow(X)
  cormat <- atanh(cormat) * sqrt(n - 3)
  d <- ncol(X)

  as.partition <- function(partition_map){
    lapply(1:max(partition_map), function(i) which(partition_map == i))
  }

  CX <- as.partition(clustering_algo(cormat, k, Xargs))
  CY <- as.partition(clustering_algo(t(cormat), k, Yargs))

  mean_cor <- function(i,j){
    mean(cormat[CX[[i]], CY[[j]]])
  }

  #Calculate mean-correlation between CX and CY.
  M <- outer(1:k, 1:k, function(x,y) mapply(mean_cor, x, y, SIMPLIFY = TRUE))

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

ircc.kmeans <- function(X, Y, k, centers.X=NULL, centers.Y=NULL) {

    clustering.kmeans <- function(data, k, args){
    if (!is.null(args)) {
      centers <- data[args,]
    } else {
      centers <- k
    }
    kmeans(data, centers, nstart = 3)$cluster
  }

  clustering_approach_to_bmd(X, Y, k, clustering.kmeans,
                                      Xargs = centers.X,
                                      Yargs = centers.Y)
}

ircc.hclust <- function(X, Y, k){
  clustering.hclust <- function(data, k, args){
    d <- dist(data)
    hc <- hclust(d)
    cutree(hc, k)
  }

  clustering_approach_to_bmd(X, Y, k, clustering.hclust)
}

ircc <- function(X, Y, nbmds, ..., method = "kmeans") {

  if (method == "kmeans") {

    res <- ircc.kmeans(X, Y, k = nbmds, ...)

  } else if (method == "hclust"){

    res <- ircc.hclust(X, Y, k = nbmds, ...)

  } else{

    stop(paste("Method", method,"not supported."))

  }

  return(res)

}
