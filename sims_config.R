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
doBRIM <- FALSE #Should we test BRIM?

saveDir <- "sims"
dataset_fname <- function(n) { file.path(saveDir, "datasets", sprintf("n=%d.RData", n)) }
results_fname <- function(n, method) { file.path(saveDir, "results", sprintf("n=%d_%s.RData", n, method)) }
plots_dir <- function(n, method){
  dirname <- file.path(saveDir, "plots", paste0("n=", n), method)
  if (!dir.exists(dirname)) dir.create(dirname, recursive = TRUE)
  return(dirname)
}
