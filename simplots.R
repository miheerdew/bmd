library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")
source("sims_config.R")

for (n in ns) {

  # Loading data and results
  load(dataset_fname(n))
  if (doBRIM) load(results_fname(n, method="brim"))
  load(results_fname(n, method="bmd"))

  #Plot the correlation matrix
  cat("computing correlation matrix")
  M <- cor(cbind(X,Y))
  cat("plotting correlation matrix")
  plot_matrix(M, filename=file.path(plots_dir(n),"corr.png"),
              width=nrow(M)/2, height=ncol(M)/2, res=400)

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
