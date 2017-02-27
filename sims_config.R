# Simulation name: corresponds to a file in "config-files" directory
# set to "" if you just want to use the default config (this file).
# Otherwise there should be a file called -
#       paste0("config-files/sims_config-",simname,".R")

simname <- "manyblocks-sparse"

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

#################################################################
#Import config-file
if(simname != ""){
  source(file.path("config-files", paste0("sims_config-", simname, ".R")))
}
