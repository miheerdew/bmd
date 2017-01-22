bX <- 50 #Size of each X block
bY <- 25 #Size of each Y block

base_rho <- seq(0.9,0,-0.1) #The intra block correlations in X.
rho_knobs <- c(0,1) #Multiplier to base_rho for the intra-block correlations

beta <- 1 #In each Y block, Y = beta * (average value of X's in block) + noise

ns <- c(100,500) #The sample sizes to run.
saveDir <- "sims"
doBRIM <- FALSE #Should we test BRIM?
