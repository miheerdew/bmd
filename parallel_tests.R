# Filematrix test

library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(rbenchmark)
library(bmdCpp)
source("makeVars.R")
source("stdize.R")

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

testnodes <- 99

res <- benchmark(Basic = {sets0 <- lapply(rep(1:testnodes, 1),
                                 function (i) initialize(n, cor(X, Y[ , i])))},
  
Par1 = {
  # 1 rep
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(bmdCpp))
  sets1 <- foreach(i = rep(1:testnodes, 1)) %dopar% {
      bmdCpp::initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
  }
  stopCluster(cl)
},
Par2 = {
  # 2 reps
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(bmdCpp))
  sets2 <- foreach(i = rep(1:testnodes, 2)) %dopar% {
    bmdCpp::initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
  }
  stopCluster(cl)
},
Par3 = {
  # 3 reps
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(bmdCpp))
  sets3 <- foreach(i = rep(1:testnodes, 3)) %dopar% {
    bmdCpp::initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
  }
  stopCluster(cl)
},
Par4 = {
  # 4 reps
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(bmdCpp))
  sets4 <- foreach(i = rep(1:testnodes, 4)) %dopar% {
    bmdCpp::initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
  }
  stopCluster(cl)
},
replications = 5)

k <- 20

res2 <- benchmark(Basic = {sets0 <- lapply(rep(1:testnodes, k),
                                 function (i) initialize(n, cor(X, Y[ , i])))},
Par = {
  # k reps
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, library(bmdCpp))
  sets1 <- foreach(i = rep(1:testnodes, k)) %dopar% {
      bmdCpp::initializeC(n, cor(X, Y[ , i]), 0.05, TRUE)
  }
  stopCluster(cl)
},
replications = 5)


