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

lambda <- 5
p <- min(lambda/nB,1) #The probability of an edge


ns <- seq(1000, 5000, by = 1000) #The sample sizes to run.
doBRIM <- FALSE #Should we test BRIM?

saveDir <- "sims"
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
