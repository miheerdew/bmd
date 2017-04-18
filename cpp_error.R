# Filematrix test

library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
source("makeVars.R")
source("stdize.R")
Rcpp::sourceCpp("bmd_helper.cpp")
no_cores <- detectCores() - 1


# Load the simulation data
load("sims-results/experiment1/90/1/sim.RData")

# Set up the data 
X <- scale(sim$X); Y <- scale(sim$Y)
X3 <- X^3; X2 <- X^2; X4ColSum <- colSums(X^4)
Y3 <- Y^3; Y2 <- Y^2; Y4ColSum <- colSums(Y^4)

dx <- ncol(X)
dy <- ncol(Y)
n  <- nrow(X)

Xindx <- 1:dx
Yindx <- (dx + 1):(dx + dy)

# Let's initialize the first 10 nodes with foreach
cl <- makeCluster(no_cores)
registerDoParallel(cl)
sets <- foreach(i = 1:10) %dopar% {
    initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
}
# Error in { : task 1 failed - "NULL value passed as symbol address"
stopCluster(cl)



# Maybe if we build the cpp functions in each node? No such luck unfortunately:
cl <- makeCluster(no_cores)
registerDoParallel(cl)
clusterEvalQ(cl, Rcpp::sourceCpp("bmd_helper.cpp"))
sets <- foreach(i = 1:10) %dopar% {
    initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
}
# Error in { : task 1 failed - "NULL value passed as symbol address"
stopCluster(cl)

