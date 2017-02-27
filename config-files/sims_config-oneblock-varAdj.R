#Number of background blocks
nBg <- 4

#Number of Bimodules
nBM <- 4

#Number of blocks per Bimodule
nB <- 1

#Block sizes
sB <- 250

# Regression coefficient
beta <- 1

# Base noise scaling
s2 <- 1

#The intra-block correlation.
rho <- 0.5


rhos <- c(rep(rho, nB*nBM), rep(0, nBg)) 
lambda <- nB
p <- min(lambda/nB,1)
eta <- ((1 - p)^nB + lambda) * sB * (1 - rho + rho * sB)
s2 <- s2 * eta
