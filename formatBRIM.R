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