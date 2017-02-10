#Block sizes
bX <- 5
bY <- 5

#Number of blocks
nX <- 20
nY <- 20

base_rho <- rep(0.9, nX) #The intra block correlations in X.
rho_knobs <- c(0,1) #Multiplier to base_rho for the intra-block correlations

beta <- 1 #In each Y block, Y = beta * (average value of X's in block) + noise

lambda <- 3
p <- min(lambda/nX,1) #The probability of an edge


ns <- c(100,500) #The sample sizes to run.
saveDir <- "sims"
doBRIM <- FALSE #Should we test BRIM?
