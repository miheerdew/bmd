library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")

# Set desired list of sample sizes to try
ns <- c(100, 500)

saveDir <- "sims"

rho_knobs <- c(0, 1)
# Do brim?
doBRIM <- FALSE

for (rcount in 1:length(rho_knobs)){
  for (ncount in 1:length(ns)) {

    fn <- paste0("rcount=", rcount, "_",
             "ncount=", ncount)

    # Loading data and results
    load(file.path(saveDir, "datasets", paste0(fn, ".RData")))
    if (doBRIM) load(file.path(saveDir, "results", paste0(fn, "_brim.RData")))
    load(file.path(saveDir, "results", paste0(fn, "_bmd.RData")))

    gen_plot_dir <- file.path("sims", "plots", paste0("rcount=",rcount), paste0("ncount=",ncount))

    if (!dir.exists(gen_plot_dir))
      dir.create(gen_plot_dir, recursive = TRUE)

    # Plotting for BMD
    sim_postanalysis(BMDresults, X, Y, run_name = paste0("BMD_",fn),
                     run_dir = file.path(gen_plot_dir, "BMD"))

    if (doBRIM) {
    # Plotting for BRIM
    sim_postanalysis(BRIMresults1, X, Y, run_name = paste0("BRIM1_", fn),
                     run_dir = file.path(gen_plot_dir, "BRIM1"), BMD = FALSE)

    sim_postanalysis(BRIMresults2, X, Y, run_name = paste0("BRIM2_", fn),
                     run_dir = file.path(gen_plot_dir, "BRIM2"), BMD = FALSE)
    }
  }
}

