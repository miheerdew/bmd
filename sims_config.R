#Block sizes
bX <- 5
bY <- 5

#Number of blocks
nX <- 20
nY <- 20

rhos <- rep(0.9, nX) #The intra block correlations in X.

beta <- 1 #In each Y block, Y = beta * (average value of X's in block) + noise

lambda <- 3
p <- min(lambda/nX,1) #The probability of an edge


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
