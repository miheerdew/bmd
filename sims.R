library(MASS)
library(Matrix)
library(lpbrim)
source("sims_config.R")

# Set intra-correlations of X's
rho_blocksX <- lapply(rhos, function (R) matrix(R, bX, bX) + diag(rep(1 - R, bX)))
mX <- nX * bX
mY <- nY * bY

#------Generating Block Covariance Matrices-----#
SigmaX <- as.matrix(bdiag(rho_blocksX))
SigmaY <- diag(mY)

for (ncount in 1:length(ns)) {

n <- ns[ncount]

set.seed(1234567)

#Random edge matrix
G <- matrix(rbinom(nY * nX, 1, p), nrow = nY, ncol = nX)
X <- mvrnorm(n, rep(0, mX), SigmaX)
Y <- mvrnorm(n, rep(0, mY), SigmaY)
Y0 <- Y

# Adding signal
for (i in 1:nY) {
  Yindices <- ((i - 1) * bY + 1):(i * bY)
  Y[ , Yindices] <- beta*as.vector(X %*% rep(G[i,],each=bX)) + Y[ , Yindices]
}

# Saving data
save(X, Y, file = dataset_fname(n))


