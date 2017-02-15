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
  
  component_list <- NULL
  bimodule_list <- NULL

  # Adding signal
  for (i in 1:nBM) {
    
    # Finding indices
    Xindices <- Yindices <- 1:(sB * nB) + (i - 1) * nB * sB
    
    # Collapsing the block effects
    effectsM <- diag(nB)[rep(1:nB, each = sB), ]
    bEffects <- beta * X[ , Xindices] %*% effectsM
    
    # Adding to the noise
    Y[ , Yindices] <- bEffects %*% G[i, , ] + Y[ , Yindices]
    
    # Finding connected components
    dimGiAdj <- nrow(G[i, , ]) + ncol(G[i, , ])
    GiAdj <- matrix(integer(dimGiAdj^2), ncol = dimGiAdj)
    crossindx1 <- 1:nrow(G[i, , ])
    crossindx2 <- (nrow(G[i, , ]) + 1):dimGiAdj
    GiAdj[crossindx1, crossindx2] <- G[i, , ]
    GiAdj[crossindx2, crossindx1] <- t(G[i, , ])
    Gi <- graph.adjacency(GiAdj, mode = "undirected")
    componentsi <- components(Gi)
    nci <- max(componentsi$membership)
    cci <- lapply(1:nci, function (c) which(componentsi$membership == c))
    
    # Adjusting indices
    cci <- lapply(cci, function (cc) {
      cc[cc > nB] <- m + Yindices[cc[cc > nB] - nB]
      ccx <- unlist(lapply(cc[cc <= nB], function (x) 1:sB + (x - 1) * sB))
      ccx <- Xindices[ccx]
      cc <- c(ccx, cc[cc > nB])
      return(cc)
      }
    )
    
    component_list <- c(component_list, cci)
    bimodule_list <- c(bimodule_list, list(c(Xindices, Yindices + m)))
    
  }
  # Saving data
  save(X, Y, component_list, bimodule_list, file = dataset_fname(n))

}
