library(MASS)
library(Matrix)
library(lpbrim)
source("sims_config.R")
source("mvrnormR.R")
source("ggcor.R")

# Set intra-correlations of X's

rho_blocksX <- lapply(rhos, function (R) matrix(R, sB, sB) + diag(rep(1 - R, sB)))
m <- (nB * nBM + nBg)*sB #The total number of X (or Y) vertices

#------Generating Block Covariance Matrices-----#
SigmaX <- as.matrix(bdiag(rho_blocksX))
SigmaY <- diag(m) * s2

for (n in ns) {

  set.seed(1234567)
  
  cat("Loop for sample size", n, "\n")

  #G[i,.,.] is the random edge matrix for the ith Bimodule
  cat("--making sub-bimodule graphs...\n")
  arrayvec <- rbinom(nBM * nB * nB * sB, 1, p)
  G <- array(arrayvec, dim = c(nBM, nB, nB * sB))
  
  # Making sure that each Y has at least 1 neighbor
  cat("--neighbor-correcting bimodule graphs...\n")
  correctY <- function (c) {
    retVec <- c
    if (sum(retVec) == 0)
      retVec[sample(length(retVec), 1)] <- 1
    return(retVec)
  }
  for (i in 1:nBM) {
    if (dim(G)[2] > 1) {
      G[i, , ] <- apply(G[i, , ], 2, correctY)
    } else {
      G[i, , ] <- apply(G[i, , , drop = FALSE], 2, correctY)
    }
  }
  
  cat("--simulating base X and Y data...\n")
  X <- mvrnormR(n, rep(0, m), SigmaX)
  Y <- mvrnormR(n, rep(0, m), SigmaY)
  
  component_list <- NULL
  bimodule_list <- NULL

  # Adding signal
  for (i in 1:nBM) {
    
    cat("--signal and components for bm", i, "...\n")
    
    # Finding indices
    Xindices <- Yindices <- 1:(sB * nB) + (i - 1) * nB * sB
    
    # Collapsing the block effects
    effectsM <- diag(nB)[rep(1:nB, each = sB), ]
    bEffects <- beta * X[ , Xindices] %*% effectsM
    
    # Adding to the noise
    Y[ , Yindices] <- bEffects %*% G[i, , ] + Y[ , Yindices]
    
    # Finding connected components
    dimGiAdj <- sum(dim(G)[2:3])
    GiAdj <- matrix(integer(dimGiAdj^2), ncol = dimGiAdj)
    crossindx1 <- 1:dim(G)[2]
    crossindx2 <- (dim(G)[2] + 1):dimGiAdj
    GiAdj[crossindx1, crossindx2] <- G[i, , ]
    GiAdj[crossindx2, crossindx1] <- t(G[i, , ])
    Gi <- igraph::graph.adjacency(GiAdj, mode = "undirected")
    componentsi <- igraph::components(Gi)
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
  
  cat("--saving data...\n")
  # Saving data
  save(X, Y, component_list, bimodule_list, file = dataset_fname(n))

  # Plotting correlation matrices
  combineData <- cbind(X, Y)
  dX <- ncol(X); dY <- ncol(Y)
  
  if (plot_full_mat) {
    cat("--calculating combined cor")
    combineCors <- cor(combineData)

    cat("--plotting full matrix...\n")
    # Order by connected components
    Xindx <- unlist(lapply(component_list, function (cc) cc[cc <= dX]))
    Yindx <- unlist(lapply(component_list, function (cc) cc[cc > dX]))
    Xbg <- setdiff(1:dX, Xindx); Xindx <- c(Xindx, Xbg)
    Ybg <- setdiff((dX + 1):(dX + dY), Yindx); Yindx <- c(Yindx, Ybg)
    combineCors <- combineCors[c(Xindx, Yindx), c(Xindx, Yindx)]
    ggcor(combineCors, file.path(plots_dir(n), "fullCors.png"), fisher = FALSE,
          title = "Full correlation matrix")
  }
  
  # Plot just connected components
  for (c in seq_along(component_list)) {
    cat("...plotting fisher values for component", c, "\n")
    cc <- component_list[[c]]
    Xindx <- cc[cc <= dX]
    Yindx <- cc[cc > dX]
    cormatc <- cor(combineData[ , c(Xindx, Yindx)])
    ggcor(cormatc, file.path(plots_dir(n), paste0("component", c, ".png")), 
          n = nrow(X), title = paste0("Component ", c, " correlations"), fisher = FALSE)
  }

}
