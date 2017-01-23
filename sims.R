library(MASS)
library(Matrix)
library(lpbrim)
source("sims_config.R")

for (rcount in 1:length(rho_knobs)) {
  
  # Set intra-correlations of X's
  rhos <- base_rho * rho_knobs[rcount]
  rho_blocksX <- lapply(rhos, function (R) matrix(R, bX, bX) + diag(rep(1 - R, bX)))
  mX <- length(rhos) * bX
  mY <- length(rhos) * bY
  
  #------Generating Block Covariance Matrices-----#
  SigmaX <- as.matrix(bdiag(rho_blocksX))
  SigmaY <- diag(mY)
  
  for (ncount in 1:length(ns)) {
    
    n <- ns[ncount]
    
    set.seed(1234567)
  
    X <- mvrnorm(n, rep(0, mX), SigmaX)
    Y <- mvrnorm(n, rep(0, mY), SigmaY)
    Y0 <- Y
    
    # Adding signal
    for (i in 1:length(rhos)) {
      Xindices <- ((i - 1) * bX + 1):(i * bX)
      Yindices <- ((i - 1) * bY + 1):(i * bY)
      Y[ , Yindices] <- as.vector(X[ , Xindices] %*% rep(beta / bX, bX)) + Y[ , Yindices]
    }
    
    # make filename
    fn <- paste0("rcount=", rcount, "_",
                 "ncount=", ncount,
                 ".RData")
    
    # Saving data
    save(X, Y, file = file.path(saveDir, "datasets", fn))
    
  }

}

