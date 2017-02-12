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
  G <- array(rbinom(nBM * nB * nB, 1, p), dim=c(nBM, nB, nB))
  X <- mvrnorm(n, rep(0, m), SigmaX)
  Y <- mvrnorm(n, rep(0, m), SigmaY)

  # Adding signal
  for (i in 1:nBM) {
    for (j in 1:nB) {
      #The indices of the jth Y block in the ith Bimodule
      Yindices <- 1:sB + sB*((j-1) + (i-1)*nB)

      #The indices corresponding to the ith Bimodule
      Xindices <- 1:(sB*nB) + (i-1)*nB*sB
      Y[ ,Yindices] <- (beta/(sB*lambda))*as.vector(X[,Xindices] %*% rep(G[i,j,],each=sB)) + Y[ , Yindices]
    }
  }

  # Saving data
  save(X, Y, file = dataset_fname(n))

}
