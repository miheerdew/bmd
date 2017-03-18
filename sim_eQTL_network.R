# Function to make eQTL-network-like data
# Argument descriptions:
# n = number of samples
# b = number of bimodules
# cmin = minimum size of bimodule half
# cmax = maximum size of bimodule half
# bgmult = number of background variables / (b * (cmax - cmin) / 2)
# betamean = mean of exponential distribution for betas
# rho = intracorrelations of X variables
# p = average edge density of bimodules
# s2 = noise variance scaling

source("mvrnormR.R")

sim_eQTL_network <- function (par_list, randomizeBeta = TRUE) {
  
  
  n = 500
  b = 10
  cmin = 50
  cmax = 50
  bgmult = 1 
  betamean = 1 
  p = 1 
  rho = 0.5
  s2 = 1
  
  avgsize <- cmin + (cmax - cmin) / 2
  
  # Setting noise variance
  mavg <- p * avgsize
  nv <- s2 * mavg * (1 - rho + rho * mavg) * betamean
  
  # Getting bimodule sizes
  if (cmin < cmax) {
    Xsizes <- sample(cmin:cmax, b, replace = TRUE)
    Ysizes <- sample(cmin:cmax, b, replace = TRUE)
  } else {
    Xsizes <- Ysizes <- rep(cmin, b)
  }
  dbg <- round(bgmult * b * avgsize)
  dx <- sum(Xsizes) + dbg; dy <- sum(Ysizes) + dbg
  
  # Initializing loop
  X <- matrix(numeric(dx * n), nrow = n); Y <- matrix(numeric(dy * n), nrow = n)
  Xindx <- 1; Yindx <- 1
  for (i in 1:b) {
    
    # Making indices
    Xindxs <- Xindx:(Xindx + Xsizes[i] - 1)
    Yindxs <- Yindx:(Yindx + Ysizes[i] - 1)
    MindxX <- 1:Xsizes[i]
    MindxY <- (Xsizes[i] + 1):(Xsizes[i] + Ysizes[i])
    
    # Generating eQTLs
    eQTLs <- matrix(rbinom(Xsizes[i] * Ysizes[i], 1, p), ncol = Ysizes[i])
    if (randomizeBeta)
      eQTLs[eQTLs == 1] <- rexp(sum(eQTLs), rate = 1 / betamean)
    
    # Making data
    SigmaX <- (1 - rho) * diag(Xsizes[i]) + 
      matrix(rho, ncol = Xsizes[i], nrow = Xsizes[i])
    Xi <- mvrnormR(n, rep(0, Xsizes[i]), SigmaX)
    noisevec <- rnorm(n * Ysizes[i], sd = sqrt(nv))
    noisemat <- matrix(noisevec, nrow = n)
    Yi <- Xi %*% eQTLs + noisemat
    
    # Putting data in, re-setting loop
    X[ , Xindxs] <- Xi; Y[ , Yindxs] <- Yi
    Xindx <- Xindx + Xsizes[i]; Yindx <- Yindx + Ysizes[i]
        
  }
  
  # Making background block
  xVars <- apply(X[ , 1:(Xindx - 1)], 2, var)
  yVars <- apply(Y[ , 1:(Yindx - 1)], 2, var)
  Xnoise <- matrix(rnorm(n * dbg, sd = sqrt(mean(xVars))), ncol = dbg)
  Ynoise <- matrix(rnorm(n * dbg, sd = sqrt(mean(yVars))), ncol = dbg)
  X[ , Xindx:dx] <- Xnoise
  Y[ , Yindx:dy] <- Ynoise
  
  return(list("X" = X, "Y" = Y))
  
}

make_param_list <- function (n = 500,
                             b = 10,
                             cmin = 100,
                             cmax = 100,
                             bgmult = 1,
                             betamean = 1,
                             p = 1,
                             rho = 0.5,
                             s2 = 1) {
  return(list("n" = n,
              "b" = b,
              "cmin" = cmin,
              "cmax" = cmax,
              "bgmult" = bgmult,
              "betamean" = betamean,
              "p" = 1,
              "rho" = rho,
              "s2" = s2))
}
    
    