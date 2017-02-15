library(MASS)
library(Matrix)
library(lpbrim)
source("sim_postanalysis.R")
source("sims_config.R")
source("ggcor.R")

for (n in ns) {

  # Loading data and results
  load(dataset_fname(n))
  if (doBRIM) load(results_fname(n, method="brim"))
  load(results_fname(n, method="bmd"))
  
  combineData <- cbind(X, Y)
  combineCors <- cor(combineData)
  
  # Order by connected components
  dX <- ncol(X); dy <- ncol(Y)
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
