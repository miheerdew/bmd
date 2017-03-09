library(MASS)
library(Matrix)
library(lpbrim)
library(igraph)
source("sim_postanalysis.R")
source("sims_config.R")
source("ggcor.R")
source("best_match.R")

# Setting up score vectors
BMD_bm <- BMD_bj <- numeric(length(ns))

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
  BMD_bm[which(ns == n)] <- bm_results['BestMatch']
  BMD_bj[which(ns == n)] <- bm_results['BackgroundMatch']
  
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
points(ns, BMD_bm, col = "red", pch = 15)
lines(ns, BMD_bm, col = "red", lty = 1)
legend("bottomleft", legend = c("BMD"), lty = 1,
       col = c("red"), pch = 15)
dev.off()


# Plotting best match scores
png(file.path(saveDir, "plots/background_jaccard.png"))
plot(0, 0, main = "Background Jaccard", xlim = range(ns) + c(-1, 1),
     ylim = c(0, 1), xlab = "n", ylab = "Background Jaccard Metric")
points(ns, BMD_bj, col = "red", pch = 15)
lines(ns, BMD_bj, col = "red", lty = 1)
legend("bottomleft", legend = c("BMD"), lty = 1,
        col = c("red"), pch = 15)
dev.off()
