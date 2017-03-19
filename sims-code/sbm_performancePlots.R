library(RColorBrewer)
library(plotrix)
library(igraph)

total_expers <- readLines("sims-results/exper-names.txt")

source("sim_eQTL_network.R")
source("best_match.R")

# Set method names:
methNames <- c("bmd")

# Set which methods to plot and their plot names
plot_meths <- c(1)
plot_names <- c("BMD")

# Set points
pchs <- c(14, 8, 3, 22, 24, 21)

# Re-get results?
getResults <- TRUE

# Plot main text?
main_text_plot <- TRUE

plot_expers <- c(1:length(total_expers))

# This should consistent throughout the experiments
nreps <- 20

# Brewing colors
colPal <- brewer.pal(9, "Set1")

# Storage for comb plots
combPlot1 <- rep(list(NULL), length(total_expers))
names(combPlot1) <- paste0("experiment", total_expers)
allscores <- array(0, dim = c(par_divs, nreps, 4))
methscores <- rep(list(allscores), length(methNames))
names(methscores) <- methNames

for (exper in plot_expers) {
  
  expString <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", expString)
  
  # Getting par values
  if (exists("axis_par_string")) {rm("axis_par_string")}
  load(file.path("sims-results/sbm-par-lists", paste0(expString, ".RData")))
  
  if (getResults) {
    
    for (p in 1:par_divs) {
      
      cat("########\n")
      cat("p=",p,"\n")
      cat("########\n")
      
      curr_dir_p <- file.path(root_dir, par_dirs[p])
      
      for (rep in 1:nreps) {
          
        cat("rep =",rep,"\n")
        curr_dir_p_rep <- file.path(curr_dir_p, rep)
        load(file.path(curr_dir_p_rep, "sim.RData"))
        dx <- ncol(sim$X); dy <- ncol(sim$Y); n <- nrow(sim$X)
                
        for (meth in methNames) {
            
          meth_fn <- file.path(curr_dir_p_rep, paste0(meth, ".RData"))
          load(meth_fn)
          
          comms <- results$communities
          full_comms <- lapply(seq_along(comms$X_sets),
                               function (i) c(comms$X_sets[[i]],                                             comms$Y_sets[[i]]))
          bg <- setdiff(1:(dx + dy), 
                        c(unlist(comms$X_sets), unlist(comms$Y_sets)))
          true_bg <- setdiff(1:(dx + dy),
                             c(unlist(sim$bm)))
          scores <- best_match_bimodule(full_comms, sim$bms,
                                        bg, true_bg)
          methscores[[meth]][p, rep, ] <- scores
        
          rm(results)
          gc()
          
        }
        
        rm(sim)
        gc()
        
      }
      
    }
    
    
      
    convertToMat <- function (resultsList, type = "mean", score = 1) {
      
      if (type == "mean") {
        stats <- unlist(lapply(resultsList, 
                               function(arr) apply(arr[ , , score], 1, mean)))
      } else {
        stats <- unlist(lapply(resultsList, 
                               function(arr) apply(arr[ , , score], 1, sd)))
      }
      
      stats <- matrix(stats, nrow = length(methNames), byrow = T)
      rownames(stats) <- methNames
      return(stats)
      
    }
      
    BM_means <- convertToMat(methscores, score = 1)
    BM1_means <- convertToMat(methscores, score = 2)
    BM2_means <- convertToMat(methscores, score = 3)
    BJ_means <- convertToMat(methscores, score = 4)
      
    BM_sds <- convertToMat(methscores, score = 1, type = "sd")
    BM1_sds <- convertToMat(methscores, score = 2, type = "sd")
    BM2_sds <- convertToMat(methscores, score = 3, type = "sd")
    BJ_sds <- convertToMat(methscores, score = 4, type = "sd")
      
    save(BM_means, BM_sds,
         BM1_means, BM1_sds,
         BM2_means, BM2_sds,
         BJ_means, BJ_sds,
         file = file.path(root_dir, "plot_results.RData"))
    
  } else {
    
    load(file.path(root_dir, "plot_results.RData"))
    
  }
  
  # Color palette
  colPal <- colPal[1:length(methNames)]
  
  source("makePerformancePlot.R")
  
  # Some plot defaults
  cex.main <- 4
  cex.lab <- 3.5
  cex.axis <- 3.5
  cex <- 3.5
  legCex <- 3
  lwd <- 3
  dotnmi <- FALSE
  
  if (exists("axis_par_string")) {
    xlab_string <- axis_par_string
  } else {
    xlab_string <- pars[axis_par]
  }
  
  paramVec <- par_settings[axis_par, ]
  
  if (paste0(expString, "_BM.png") %in% 
      list.files(file.path("sims-results", expString))) {
    file.remove(paste0(expString, "_BM.png"))
  }

  if (main_text_plot) {
    main_str <- main_text
  } else {
    main_str <- ""
  }
  
  
  BM_means <- BM_means[plot_meths, , drop = FALSE]
  BM_sds <- BM_sds[plot_meths, , drop = FALSE]
  rownames(BM_means) <- plot_names
  rownames(BM_sds) <- plot_names
  
  BM1_means <- BM1_means[plot_meths, , drop = FALSE]
  BM1_sds <- BM1_sds[plot_meths, , drop = FALSE]
  rownames(BM1_means) <- plot_names
  rownames(BM1_sds) <- plot_names
  
  BM2_means <- BM2_means[plot_meths, , drop = FALSE]
  BM2_sds <- BM2_sds[plot_meths, , drop = FALSE]
  rownames(BM2_means) <- plot_names
  rownames(BM2_sds) <- plot_names
  
  BJ_means <- BJ_means[plot_meths, , drop = FALSE]
  BJ_sds <- BJ_sds[plot_meths, , drop = FALSE]
  rownames(BJ_means) <- plot_names
  rownames(BJ_sds) <- plot_names
  
  plotfn = file.path("sims-results", paste0(expString, ".png"))
  
  png(plotfn, width = 1000, height = 1000)
  par(mfrow = c(2, 2), oma = rep(0, 4),
      mar = c(11, 11, 6, 4),
      mgp = c(6, 2, 0))
  
  # BM
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, 
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM_means, tnmi = dotnmi,
                                 sdMat = BM_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Jaccard",
                                 legPos = "bottomright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
  )
    
  # BM1
  
  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM1_means, tnmi = dotnmi,
                                 sdMat = BM1_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Res.prop",
                                 legPos = "bottomright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # BM2

  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BM2_means, tnmi = dotnmi,
                                 sdMat = BM2_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Best Match Truth.prop",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  # BJ

  suppressWarnings(
    
    dummy <- makePerformancePlot(plotFile = FALSE, doLegend = FALSE,
                                 xvals = paramVec, yRange = c(0, 1),
                                 meanMat = BJ_means, tnmi = dotnmi,
                                 sdMat = BJ_sds,
                                 xRange = c(paramVec[1], paramVec[length(paramVec)]),
                                 main = main_str,
                                 xlab = xlab_string,
                                 ylab = "Background Jaccard",
                                 legPos = "topright",
                                 legCex = legCex,
                                 lwd = lwd,
                                 cex = cex, pchs = pchs,
                                 cex.main = cex.main,
                                 cex.lab = cex.lab,
                                 cex.axis = cex.axis)
    
  )
  
  dev.off()
  
}