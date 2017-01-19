library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")

# Set desired list of sample sizes to try
ns <- c(100, 500)

saveDir <- "sims"

# Do brim?
doBRIM <- FALSE

for (ncount in 1:length(ns)) {
  
  # Loading data and results
  load(file.path(saveDir, "datasets", paste0(ncount, ".RData")))
  load(file.path(saveDir, "results", paste0("ncount", ncount, "_brim.RData")))
  load(file.path(saveDir, "results", paste0("ncount", ncount, "_bmd.RData")))
  
  gen_plot_dir <- file.path("sims", "plots", ncount)
  
  if (!dir.exists(gen_plot_dir))
    dir.create(gen_plot_dir, recursive = TRUE)
  
  # Plotting for BMD
  sim_postanalysis(BMDresults, X, Y, run_name = paste0("BMD_ncount=", ncount),
                   run_dir = file.path(gen_plot_dir, "BMD"))
  
  if (doBRIM) {
  # Plotting for BRIM
  sim_postanalysis(BRIMresults1, X, Y, run_name = paste0("BRIM1_ncount=", ncount),
                   run_dir = file.path(gen_plot_dir, "BRIM1"), BMD = FALSE)
  
  sim_postanalysis(BRIMresults2, X, Y, run_name = paste0("BRIM2_ncount=", ncount),
                   run_dir = file.path(gen_plot_dir, "BRIM2"), BMD = FALSE)
  }
  
}


