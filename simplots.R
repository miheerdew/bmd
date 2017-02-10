library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")
source("sims_config.R")

for (rcount in 1:length(rho_knobs)){
  for (ncount in 1:length(ns)) {
    n <- ns[ncount]

    # Loading data and results
    load(dataset_fname(n))
    if (doBRIM) load(results_fname(n, method="brim"))
    load(results_fname(n, method="bmd"))

    fn = sprintf("n=%d", n)

    # Plotting for BMD
    sim_postanalysis(BMDresults, X, Y, run_name = paste0("BMD_",fn),
                     run_dir = plots_dir(n, method="BMD"))

    if (doBRIM) {
    # Plotting for BRIM
    sim_postanalysis(BRIMresults1, X, Y, run_name = paste0("BRIM1_", fn),
                     run_dir = plots_dir(n, method="BRIM1"), BMD = FALSE)

    sim_postanalysis(BRIMresults2, X, Y, run_name = paste0("BRIM2_", fn),
                     run_dir = plots_dir(n, method="BRIM2"), BMD = FALSE)
    }
  }
}

