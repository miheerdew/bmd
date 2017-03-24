library(lpbrim)

formatBRIM <- function (mod, zeroDegs, dx) {
  
  ncomms <- ncol(mod[[1]])
  p <- length(zeroDegs)
  dxp <- sum(which(!zeroDegs) <= dx)
  fullmods <- lapply(1:ncomms, function (c) rep(which(!zeroDegs), 2)[as.logical(mod[[1]][  , c])])
  side1mods <- unlist(lapply(fullmods, function (L) L[1] <= dxp))
  
  X_sets1 <- lapply(fullmods[side1mods], function (L) L[L <= dx])
  X_sets2 <- lapply(fullmods[!side1mods], function (L) L[L <= dx])
  Y_sets1 <- lapply(fullmods[side1mods], function (L) L[L > dx])
  Y_sets2 <- lapply(fullmods[!side1mods], function (L) L[L > dx])
  
  Xsizes1 <- unlist(lapply(X_sets1, length))
  Ysizes1 <- unlist(lapply(Y_sets1, length))
  singletons1 <- Ysizes1 + Xsizes1 <= 2
  X_sets1 <- X_sets1[!singletons1]
  Y_sets1 <- Y_sets1[!singletons1]
  
  Xsizes2 <- unlist(lapply(X_sets2, length))
  Ysizes2 <- unlist(lapply(Y_sets2, length))
  singletons2 <- Ysizes2 + Xsizes2 <= 2
  X_sets2 <- X_sets2[!singletons2]
  Y_sets2 <- Y_sets2[!singletons2]
  
  
  communities1 <- list("X_sets" = X_sets1,
                       "Y_sets" = Y_sets1)
  
  communities2 <- list("X_sets" = X_sets2,
                       "Y_sets" = Y_sets2)
  return(list(communities1, communities2))
  
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
  largeM <- matrix(integer(p^2), ncol = p)
  largeM[1:p1, (p1 + 1):p] <- crossSigMat
  largeM[(p1 + 1):p, 1:p1] <- t(crossSigMat)
  
  # Making sure no rowSums are zero
  zeroDegs <- rowSums(largeM) == 0
  largeM_red <- largeM[!zeroDegs, !zeroDegs]
  if (sum(zeroDegs) == length(zeroDegs)) {
    BRIMresults1 <- BRIMresults2 <- 
      list("communities" = list("X_sets" = NULL, "Y_sets" = NULL))
  } else {
    BRIMmod <- findModules(largeM_red, iter = 5, sparse = FALSE)
    BRIMformat <- formatBRIM(BRIMmod, zeroDegs, ncol(X))
    BRIMresults1 <- list("communities" = BRIMformat[[1]])
    BRIMresults2 <- list("communities" = BRIMformat[[2]])
  }
  
  return(list("BRIMresults1" = BRIMresults1,
              "BRIMresults2" = BRIMresults2))
  
}