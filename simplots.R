library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")
source("sims_config.R")
source("ggcor.R")
source("best_match.R")

# Setting up score vectors
BMD_bm <- BMD_bj <- rep(0, length(ns))

plot_full_mat <- FALSE

for (n in ns) {

  # Loading data and results
  load(dataset_fname(n))
  if (doBRIM) load(results_fname(n, method="brim"))
  load(results_fname(n, method="bmd"))
  
  combineData <- cbind(X, Y)
  combineCors <- cor(combineData)
  dX <- ncol(X); dY <- ncol(Y)
  
  if (plot_full_mat) {
    # Order by connected components
    Xindx <- unlist(lapply(component_list, function (cc) cc[cc <= dX]))
    Yindx <- unlist(lapply(component_list, function (cc) cc[cc > dX]))
    Xbg <- setdiff(1:dX, Xindx); Xindx <- c(Xindx, Xbg)
    Ybg <- setdiff((dX + 1):(dX + dY), Yindx); Yindx <- c(Yindx, Ybg)
    combineCors <- combineCors[c(Xindx, Yindx), c(Xindx, Yindx)]
    cat("...plotting full correlation matrix\n")
    ggcor(combineCors, file.path(plots_dir(n), "fullCors.png"), fisher = FALSE,
          title = "Full correlation matrix")
    cat("...plotting full fisher value matrix\n")
    ggcor(combineCors, file.path(plots_dir(n), "fullFish.png"), n = nrow(X),
          title = "Full fisher transform matrix")
  }
  
  # Plot just connected components
  for (c in seq_along(component_list)) {
    cat("...plotting fisher values for component", c, "\n")
    cc <- component_list[[c]]
    Xindx <- cc[cc <= dX]
    Yindx <- cc[cc > dX]
    cormatc <- cor(combineData[ , c(Xindx, Yindx)])
    ggcor(cormatc, file.path(plots_dir(n), paste0("component", c, ".png")), 
          n = nrow(X), title = paste0("Component ", c, " fisher values"))
  }

  fn = sprintf("n=%d", n)
  
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
  BMDresults2 <- list(BMDresults_comms, 
                      c(BMDresults$background$X_bg,
                        BMDresults$background$Y_bg))
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
