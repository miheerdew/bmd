library(MASS)
library(Matrix)
library(lpbrim)
library(igraph)
source("sim_postanalysis.R")
source("sims_config.R")
source("ggcor.R")
source("best_match.R")

# Setting up score vectors
BMD_bm <- BMD_bj <- matrix(0, 2, length(ns))

for (n in ns) {
  
  cat("plotting for n =", n, "\n")

  fn = sprintf("n=%d", n)
  
  # Loading data and results
  load(dataset_fname(n))

  # Plotting for BMD
  cat("----doing BMD\n")
  load(results_fname(n, method="bmd"))
  bm_results <- sim_postanalysis(BMDresults, X, Y, run_name = paste0("bmd_",fn),
                   run_dir = plots_dir(n, method="bmd"))
  BMD_bm[1, which(ns == n)] <- bm_results['BestMatch']
  BMD_bj[1, which(ns == n)] <- bm_results['BackgroundMatch']
  
  rm(bm_results, BMDtime, BMDresults)
  gc()
  
  # Plotting for BMD-kb
  cat("----doing BMD-kb\n")
  load(results_fname(n, method="bmd-kb"))
  bm_results <- sim_postanalysis(BMDresults, X, Y, run_name = paste0("bmd_",fn),
                   run_dir = plots_dir(n, method="bmd-kb"))
  BMD_bm[2, which(ns == n)] <- bm_results['BestMatch']
  BMD_bj[2, which(ns == n)] <- bm_results['BackgroundMatch']
  
  rm(bm_results, BMDtime, BMDresults)
  gc()

  if (doBRIM) {
    if (doBRIM) load(results_fname(n, method="brim"))
    # Plotting for BRIM
    sim_postanalysis(BRIMresults1, X, Y, run_name = paste0("BRIM1_", fn),
                     run_dir = plots_dir(n, method="BRIM1"), BMD = FALSE)

    sim_postanalysis(BRIMresults2, X, Y, run_name = paste0("BRIM2_", fn),
                     run_dir = plots_dir(n, method="BRIM2"), BMD = FALSE)
  }
  
  
}

# Plotting best match scores
png(file.path(saveDir, "plots/best_match.png"))
plot(0, 0, main = "Best Match", xlim = range(ns) + c(-1, 1),
     ylim = c(0, 1), xlab = "n", ylab = "Best Match Metric")
points(ns, BMD_bm[1, ], col = "red", pch = 15)
lines(ns, BMD_bm[1, ], col = "red", lty = 1)
points(ns, BMD_bm[2, ], col = "green", pch = 16)
lines(ns, BMD_bm[2, ], col = "green", lty = 2)
legend("bottomleft", legend = c("BMD", "BMD-kb"), lty = 1:2,
       col = c("red", "green"), pch = 15:16)
dev.off()


# Plotting best match scores
png(file.path(saveDir, "plots/background_jaccard.png"))
plot(0, 0, main = "Background Jaccard", xlim = range(ns) + c(-1, 1),
     ylim = c(0, 1), xlab = "n", ylab = "Background Jaccard Metric")
points(ns, BMD_bj[1, ], col = "red", pch = 15)
lines(ns, BMD_bj[1, ], col = "red", lty = 1)
points(ns, BMD_bj[2, ], col = "green", pch = 16)
lines(ns, BMD_bj[2, ], col = "green", lty = 2)
legend("bottomleft", legend = c("BMD", "BMD-kb"), lty = 1:2,
       col = c("red", "green"), pch = 15:16)
dev.off()
