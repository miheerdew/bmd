# Name of sim
simname <- "oneblock-varAdj"

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


ns <- c(100, 500, seq(1000, 5000, by = 1000)) #The sample sizes to run.
doBRIM <- FALSE #Should we test BRIM?

saveDir <- file.path("sims", simname)
dataset_fname <- function(n) {
  if (!dir.exists(file.path(saveDir, "datasets")))
    dir.create(file.path(saveDir, "datasets"), recursive = TRUE)
  file.path(saveDir, "datasets", sprintf("n=%d.RData", n))
}
results_fname <- function(n, method) {
  if (!dir.exists(file.path(saveDir, "results")))
    dir.create(file.path(saveDir, "results"), recursive = TRUE)
  file.path(saveDir, "results", sprintf("n=%d_%s.RData", n, method))
}
plots_dir <- function(n, method=""){
  dirname <- file.path(saveDir, "plots", paste0("n=", n), method)
  if (!dir.exists(dirname)) dir.create(dirname, recursive = TRUE)
  return(dirname)
}
