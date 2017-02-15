library(MASS)
library(Matrix)
library(lpbrim)
source("sims_config.R")

# Set intra-correlations of X's

rho_blocksX <- lapply(rhos, function (R) matrix(R, sB, sB) + diag(rep(1 - R, sB)))
m <- (nB * nBM + nBg)*sB #The total number of X (or Y) vertices

#------Generating Block Covariance Matrices-----#
SigmaX <- as.matrix(bdiag(rho_blocksX))
SigmaY <- diag(m)

for (n in ns) {

  set.seed(1234567)

  #G[i,.,.] is the random edge matrix for the ith Bimodule
  G <- array(rbinom(nBM * nB * nB * sB, 1, p), dim = c(nBM, nB, nB * sB))
  
  # Making sure that each Y has at least 1 neighbor
  correctY <- function (c) {
    retVec <- c
    if (sum(retVec) == 0)
      retVec[sample(length(retVec), 1)] <- 1
    return(retVec)
  }
  for (i in 1:nBM) {
    G[i, , ] <- apply(G[i, , ], 2, correctY)
  }
  X <- mvrnorm(n, rep(0, m), SigmaX)
  Y <- mvrnorm(n, rep(0, m), SigmaY)

  # Adding signal
  for (i in 1:nBM) {
    
    # Finding indices
    Xindices <- Yindices <- 1:(sB * nB) + (i - 1) * nB * sB
    
    # Collapsing the block effects
    effectsM <- diag(nB)[rep(1:nB, each = sB), ]
    bEffects <- beta * X[ , Xindices] %*% effectsM
    
    # Adding to the noise
    Y[ , Yindices] <- bEffects %*% G[i, , ]
    
  }
  # Saving data
  save(X, Y, file = dataset_fname(n))

}
