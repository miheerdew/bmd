library(MASS)
library(Matrix)
library(lpbrim)
library(igraph)
source("sim_postanalysis.R")
source("sims_config.R")
source("ggcor.R")
source("best_match.R")

# Setting up score vectors
BMD_bm <- BMD_bj <- rep(0, length(ns))

for (n in ns[-1]) {

  # Loading data and results
  load(dataset_fname(n))
  if (doBRIM) load(results_fname(n, method="brim"))
  load(results_fname(n, method="bmd"))
  
  fn = sprintf("n=%d", n)
  
  dX <- ncol(X)
  dY <- ncol(Y)
  
  # Preparing component list for scoring
  component_list2 <- list(component_list,
                          setdiff(1:(dX + dY), unlist(component_list)))
  
  # Preparing BMD results for scoring
  BMDresults_comms <- lapply(seq_along(BMDresults$communities[[1]]), 
                             function (c) {
                               c(BMDresults$communities$X_sets[[c]],
                                 BMDresults$communities$Y_sets[[c]])
                             }
  )
  BMDresults2 <- c(list(BMDresults_comms), 
                   list(c(BMDresults$background$X_bg,
                          BMDresults$background$Y_bg)))
  BMD_bestmatch <- best_match_bimodule(component_list2, BMDresults2)
  BMD_bm[which(ns == n)] <- BMD_bestmatch['BestMatch']
  BMD_bj[which(ns == n)] <- BMD_bestmatch['BackgroundMatch']

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

# Plotting best match scores
png("sims/plots/best_match.png")
plot(ns, BMD_bm, main = "Best Match")
dev.off()


# Plotting background jaccard scores
png("sims/plots/background_jaccard.png")
plot(ns, BMD_bj, main = "Background Jaccard")
dev.off()
