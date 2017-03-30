library(lpbrim)

formatBRIM <- function (mod, zeroDegsR, zeroDegsC, dx) {
  
  ncomms <- ncol(mod[[1]])
  dx <- length(zeroDegsR)
  dy <- length(zeroDegsC)
  goodnodes <- c(which(!zeroDegsR), (c((dx + 1):(dx + dy)))[which(!zeroDegsC)])
  fullmods <- lapply(1:ncomms, function (c) goodnodes[as.logical(mod[[1]][  , c])])
  fullmods <- fullmods[unlist(lapply(fullmods, length)) > 2]
  X_sets <- lapply(fullmods, function (C) C[C <= dx])
  Y_sets <- lapply(fullmods, function (C) C[C > dx])
  retlist <- list("communities" = list("X_sets" = X_sets,
                                       "Y_sets" = Y_sets),
                  "background" = setdiff(1:(dx + dy), union(unlist(X_sets),
                                                            unlist(Y_sets))))
  return(retlist)
  
}

run_brim <- function(X, Y, alpha = 0.05) {
  
  # Making large correlation matrix
  n <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Y)
  p <- p1 + p2
  crossCors <- cor(X, Y)
  crossTstats <- crossCors * sqrt(n - 2) / sqrt(1 - crossCors^2)
  crossPvals <- as.vector(pt(crossTstats, n - 2, lower.tail = FALSE))
  crossPvals_adj <- crossPvals * p1 * p2 / rank(crossPvals)
  thres <- max(crossPvals[which(crossPvals_adj <= alpha)])
  crossSigs <- as.integer(crossPvals_adj <= thres)
  crossSigMat <- matrix(crossSigs, ncol = p2)
  
  # Making sure no r/c are zero
  zeroDegsR <- rowSums(crossSigMat) == 0
  zeroDegsC <- colSums(crossSigMat) == 0
  CSM_red <- crossSigMat[!zeroDegsR, !zeroDegsC]
  if (sum(zeroDegs) == length(zeroDegs)) {
    BRIMresults <- list("communities" = list("X_sets" = NULL, "Y_sets" = NULL),
                        "background" = 1:(p1 + p2))
  } else {
    BRIMmod <- findModules(CSM_red, iter = 5, sparse = FALSE)
    BRIMresults <- formatBRIM(BRIMmod, zeroDegsR, zeroDegsC, ncol(X))
  }
  
  return(BRIMresults)
  
}