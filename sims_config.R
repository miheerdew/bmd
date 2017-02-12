#Number of background blocks
nBg <- 200

#Number of Bimodules
nBM <- 4

#Number of blocks per Bimodule
nB <- 50

#Block sizes
sB <- 5

rhos <- c(rep(0.5, nB*nBM), rep(0, nBg)) #The intra-block correlation.

beta <- 1

lambda <- 3
p <- min(lambda/nB,1) #The probability of an edge


ns <- c(100,500) #The sample sizes to run.
saveDir <- "sims"
doBRIM <- FALSE #Should we test BRIM?
