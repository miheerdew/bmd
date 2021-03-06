library(Matrix)
source("plotLargeMat.R")

plot_matrix <- function(m, filename, width=600, height=600,
                        useRaster=FALSE, ... ) {
  #Plots the square matrix to a png file.
  png(filename, width=width, height=height, ...)
  image(m[,ncol(m):1], useRaster=useRaster)
  dev.off()
}

sim_postanalysis <- function (results, X, Y, run_name, run_dir,
                              X_name = "X", Y_name = "Y", BMD = TRUE, do_update_info = TRUE,
                              GTEx = FALSE, snpsloc = NULL, geneloc = NULL) {

  #-----------------------------------------------------------------------

  if (!dir.exists(run_dir)) {
    dir.create(run_dir)
  }

  cat("calculating full cross-correlation matrix...\n")
  full_cross_cor <- cor(X,Y)

  cat("calculating basic stats...\n")

  # Calculating basic stats
  ncomms <- length(results$communities$X_sets)
  comm_dfs <- list(rep(list(NULL), ncomms))
  Y_ids <- (ncol(X) + 1):(ncol(X) + ncol(Y))

  # Getting scaled X and Y
  X_scl <- scale(X)
  Y_scl <- scale(Y)

  Xsizes <- integer(ncomms)
  Ysizes <- integer(ncomms)
  
  if (do_update_info) {
    
    nonNullUpdateIndxs <- which(!unlist(lapply(results$update_info, length)) == 0)
    nonNullIndxs <- results$nonNullIndxs
    commIndxs <- results$finalIndxs
    
    # Creating folder for update info plots
    uiplotdir <- file.path(run_dir, "update_info")
    if (!dir.exists(uiplotdir))
      dir.create(uiplotdir)
    
    for (ui in nonNullUpdateIndxs) {
      
      updateInfoi <- results$update_info[ui]
      
      uilength <- length(updateInfoi[[1]]$consec_jaccards)
      
      plottitle <- paste0("ui ", ui, ",")
      fn <- paste0("ui", ui)
      if (ui %in% commIndxs) {
        plottitle <- paste(plottitle, "found community")
        fn <- paste0(fn, "fc")
      }
      if (ui %in% setdiff(nonNullIndxs, commIndxs)) {
        plottitle <- paste(plottitle, "overlapping community")
        fn <- paste0(fn, "oc")
      }
      if (ui %in% setdiff(nonNullUpdateIndxs, nonNullIndxs)) {
        plottitle <- paste(plottitle, "broken cycle or null community")
        fn <- paste0(fn, "bc")
      }
      
      png(file.path(uiplotdir, paste0(fn, ".png")))
      par(mfrow=c(1,2))

      #Plot the consecutive jaccards for the extraction
      plot(0, 0, xlim = c(0, uilength + 1), ylim = c(0, 1), col = "white",
           main = plottitle, xlab = "update no.", ylab = "consec. jaccard")
      lines(1:uilength, updateInfoi[[1]]$consec_jaccards)
      points(1:uilength, updateInfoi[[1]]$consec_jaccards,
             col = as.numeric(updateInfoi[[1]]$found_cycle) + 1,
             pch = as.numeric(updateInfoi[[1]]$found_break) + 16)
      legend("topright", legend = c("found cycle", "found break"),
             col = c(2, 1), lty = c(1, NA), pch = c(NA, 17), lwd = c(4, NA))

      #Plot the consecutive module sizes for the extraction
      consec_size_pairs <- updateInfoi[[1]]$consec_sizes
      y <- total_sizes <- sapply(consec_size_pairs, sum) #total module size
      plot(0, 0, xlim = c(0, uilength + 1), ylim = range(y) + c(-1,1), col = "white",
            main = "module sizes", xlab = "update no.", ylab = "Module size")
      lines(0:uilength, y)
      points(0:uilength, y,
             col = as.numeric(updateInfoi[[1]]$found_cycle) + 1,
             pch = as.numeric(updateInfoi[[1]]$found_break) + 16)
      dev.off()
    }

  }

  for (c in 1:ncomms) {

    cat(c, "/", ncomms, "\n")

    nX <- length(results$communities$X_sets[[c]])
    nY <- length(results$communities$Y_sets[[c]])
    Xsizes[c] <- nX
    Ysizes[c] <- nY
    Y_mch <- match(results$communities$Y_sets[[c]], Y_ids)
    X_mch <- results$communities$X_sets[[c]]

    comm_df <- data.frame("X_indx" = rep(X_mch, each = nY),
                          "Y_indx" = rep(Y_mch, nX))

    corrs <- full_cross_cor[X_mch, Y_mch]
    comm_df$corrs <- as.vector(t(corrs))

    if (GTEx) {
      X_snps <- match(S[X_mch], snpsloc$ID)
      Y_gene <- match(G[Y_mch], geneloc$Id)
      X_pos  <- snpsloc$POS[X_snps]
      Y_pos  <- geneloc$TSS[Y_gene]
      comm_df$X_pos <- rep(X_pos, each = nY)
      comm_df$Y_pos <- rep(Y_pos, nX)
      comm_df$snp <- rep(S[X_mch], each = nY)
      comm_df$gene <- rep(G[Y_mch], nX)
    }

    comm_dfs[[c]] <- comm_df

  }

  commsizes <- unlist(lapply(comm_dfs, nrow))
  npairs <- sum(commsizes)
  n <- nrow(X)


  #-------------------------------------------------------------
  # COMPUTE: Community size vs -log10 pvalue distribution
  #-------------------------------------------------------------
  cat("calculating t-stats and p-values...\n")
  tstats <- unlist(lapply(comm_dfs, function (df) df$corrs * sqrt(n - 2) / sqrt(1 - df$corrs^2)))
  log10pvals <- -log10(pt(tstats, n - 2, lower.tail = FALSE))
  sizes <- unlist(lapply(commsizes, function (m) rep(m, m)))

  #-------------------------------------------------------------
  # PLOT: Community size vs -log10 pvalue distribution
  #-------------------------------------------------------------
  cat("plotting comm size vs. p-values...\n")
  fn <- "sizeVsPvals"

  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)

  boxplot(log10pvals~round(log10(sizes), 2), outline = FALSE,
          xlab = "# suggested pairs (log10 scale, numbers not aligned to axis)",
          ylab = "boxplot of t p-values (-log10 scale)",
          main = paste0("Sugg. pairs vs. -log10 p-values, ", run_name))

  dev.off()

  #-------------------------------------------------------------

  if (GTEx) {

    #-------------------------------------------------------------
    # COMPUTE: Community size vs trans distance distribution
    #-------------------------------------------------------------
    cat("calculating eQTL distances...\n")
    dists <- unlist(lapply(comm_dfs, function (df) df$X_pos - df$Y_pos))

    #-------------------------------------------------------------
    # PLOT: Community size vs trans distance distribution
    #-------------------------------------------------------------
    cat("plotting eQTL distances...\n")
    fn <- "sizeVsDist"

    png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)

    boxplot(abs(dists)~round(log10(sizes), 2), outline = FALSE,
            xlab = "# suggested eQTLs (log10 scale, numbers not aligned to axis)",
            ylab = "boxplot of genome distances (log10 scale)",
            main = paste0("Sugg. eQTLs vs. genome distance, ", run_name))

    dev.off()

    fn_string <- c(fn_string, fn)

    #-------------------------------------------------------------


    #-------------------------------------------------------------
    # PLOT: Median community trans distance vs -log10 pvalue distribution
    #-------------------------------------------------------------
    cat("calculating median absolute eQTL distances...\n")
    med_abs_dists <- unlist(lapply(comm_dfs, function (df) rep(median(abs(df$X_pos - df$Y_pos)), nrow(df))))

    #-------------------------------------------------------------
    # PLOT: Median community trans distance vs -log10 pvalue distribution
    #-------------------------------------------------------------
    cat("plotting median absolute eQTL distances...\n")
    fn <- "med_distVsPvals"

    png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)

    boxplot(log10pvals~round(log10(med_abs_dists), 2), outline = FALSE,
            xlab = "median of genome distances (log10 scale, numbers not aligned to axis)",
            ylab = "boxplot of p-values (-log10 scale)",
            main = paste0("Median genome distance vs. -log10 p-values, ", run_name))

    dev.off()

    #-------------------------------------------------------------

  }


  #-------------------------------------------------------------
  # PLOT: SNP count vs gene count
  #-------------------------------------------------------------
  cat("plotting Xcount vs Ycount...\n")
  fn <- "XcountVsYcount"

  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)

  plot(Xsizes, Ysizes, col = "deepskyblue4", pch = 18, cex = 2,
       xlab = paste0("# of ", X_name),
       ylab = paste0("# of ", Y_name),
       main = paste0("# ", X_name, " vs. # ", Y_name, ", ", run_name))

  dev.off()


  #-------------------------------------------------------------
  # COMPUTE: correlation distributions
  #-------------------------------------------------------------
  cat("computing within-community correlations...\n")
  comm_corrs <- unlist(lapply(comm_dfs, function (df) df$corrs))
  cat("computing all correlations...\n")
  all_corrs <- as.vector(full_cross_cor)

  cat("formatting for correlation density plot...\n")
  corrs0 <- density(all_corrs)
  corrs1 <- density(comm_corrs)
  plot_max <- max(c(corrs0$y, corrs1$y))
  plotlim_x <- max(abs(min(c(corrs0$x, corrs1$x))), max(c(corrs0$x, corrs1$x)))
  cols <- c("black", "deeppink")

  #-------------------------------------------------------------
  # PLOT: correlation distributions
  #-------------------------------------------------------------
  cat("plotting correlation densities...\n")
  fn <- "corr_dists"

  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)
  plot(0, 0, col = "white", xlim = c(-plotlim_x, plotlim_x), ylim = c(0, plot_max),
       main = "Correlation densities",
       xlab = "Pearson correlation",
       ylab = "Density")
  lines(corrs0, lwd = 3, col = cols[1])
  lines(corrs1, lwd = 3, col = cols[2])
  legend("topleft", c("Population", "Modules"), col = cols, lwd = 3)
  dev.off()

  #-------------------------------------------------------------
  # COMPUTE: commsizes vs median corr
  #-------------------------------------------------------------
  meancorrs <- unlist(lapply(comm_dfs, function (L) mean(L$corrs)))

  #-------------------------------------------------------------
  # PLOT: commsizes vs mean corr
  #-------------------------------------------------------------
  cat("plotting commsizes vs mean corr...\n")

  fn <- "Xsizes_vs_mcorr"
  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)
  plot(Xsizes, meancorrs, col = "deepskyblue4", pch = 18, cex = 2,
       main = "Mean correlation vs. # X",
       xlab = paste0("# of ", X_name),
       ylab = "Mean correlation")
  dev.off()

  fn <- "Ysizes_vs_mcorr"
  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)
  plot(Ysizes, meancorrs, col = "deepskyblue4", pch = 18, cex = 2,
       main = "Mean correlation vs. # Y",
       xlab = paste0("# of ", Y_name),
       ylab = "Mean correlation")
  dev.off()


  #-------------------------------------------------------------
  # COMPUTE: overlap
  #-------------------------------------------------------------
  cat("computing overlap...\n")
  X_overlap <- matrix(0, ncomms, ncomms)
  Y_overlap <- X_overlap
  for (i in 1:ncomms) {

    for (j in 1:ncomms) {

      X_i_in_j <- length(intersect(results$communities$X_sets[[i]], results$communities$X_sets[[j]])) / length(results$communities$X_sets[[i]])
      Y_i_in_j <- length(intersect(results$communities$Y_sets[[i]], results$communities$Y_sets[[j]])) / length(results$communities$Y_sets[[i]])
      X_overlap[i, j] <- X_i_in_j
      Y_overlap[i, j] <- Y_i_in_j

    }
  }


  #-------------------------------------------------------------
  # PLOT: overlap
  #-------------------------------------------------------------
  cat("plotting overlapp...\n")
  fn <- "comm_overlap"

  png(file.path(run_dir, paste0(run_name, "_", fn, "_", X_name, ".png")), width = 600, height = 600)
  image(as.matrix(X_overlap[ncol(X_overlap):1, ]))
  dev.off()

  png(file.path(run_dir, paste0(run_name, "_", fn, "_", Y_name, ".png")), width = 600, height = 600)
  image(as.matrix(Y_overlap[ncol(Y_overlap):1, ]))
  dev.off()

  #-------------------------------------------------------------
  # COMPUTE: commsizes
  #-------------------------------------------------------------
  commsizes_sum <- Xsizes + Ysizes
  commcorrs_sum <- unlist(lapply(comm_dfs, function (L) sum(L$corrs)))
  commsizes_ord <- order(commcorrs_sum, decreasing = TRUE)

  #-------------------------------------------------------------
  # PLOT: commsizes
  #-------------------------------------------------------------
  cat("plotting commsizes...\n")

  fn <- "commsizes"
  png(file.path(run_dir, paste0(run_name, "_", fn, ".png")), width = 600, height = 600)
  plot(commsizes_sum[commsizes_ord], col = "deepskyblue4", pch = 18, cex = 2,
       main = "Module sizes (sum)",
       xlab = paste0("Comm #"),
       ylab = paste0("# ", X_name, " + # ", Y_name))
  dev.off()

  #-------------------------------------------------------------
  # COMPUTE: CORMAT
  #-------------------------------------------------------------
  crossCors <- full_cross_cor
  Yindxs <- nrow(crossCors) + 1:ncol(crossCors)
  X_sets <- results$communities$X_sets[commsizes_ord]
  Y_sets <- lapply(results$communities$Y_sets, function (C) match(C, Yindxs))[commsizes_ord]

  #-------------------------------------------------------------
  # PLOT: CORMAT
  #-------------------------------------------------------------
  fn <- "cormat_clusters"
  png(file.path(run_dir, paste0(run_name, "_", fn, "_", X_name, ".png")), width = 600, height = 600)
  plotLargeMat(crossCors, X_sets, Y_sets, maxPlot = 1e6)
  dev.off()

  #-------------------------------------------------------------
  # COMPUTE: ACCEPTANCE CURVES
  #-------------------------------------------------------------

  if (BMD) {
    #-------------------------------------------------------------
    # PLOT: ACCEPTANCE CURVES
    #-------------------------------------------------------------
    fn <- "acceptance_curves"
    if (!dir.exists(file.path(run_dir, "acceptance_curves"))) {
      dir.create(file.path(run_dir, "acceptance_curves"))
    } else {
      acc_dir_list <- list.files(file.path(run_dir, "acceptance_curves"), full.names = TRUE)
      sapply(acc_dir_list, file.remove)
    }
    for (c in 1:ncomms) {

      XL <- length(results$communities$X_sets[[c]])
      YL <- length(results$communities$Y_sets[[c]])

      png(file.path(run_dir, "acceptance_curves",
                    paste0(run_name, "_", fn, "_", X_name, "_comm", c, ".png")),
          width = 600, height = 600)
      plot(results$commzs$X_zs[[c]]$zs,
           main = paste0(X_name, " module ", c, ", acceptance curve"),
           xlab = "Index", ylab = "z-statistic")
      points(results$commzs$X_zs[[c]]$zs[1:XL], col = "red")
      dev.off()

      png(file.path(run_dir, "acceptance_curves",
                    paste0(run_name, "_", fn, "_", Y_name, "_comm", c, ".png")),
          width = 600, height = 600)
      plot(results$commzs$Y_zs[[c]]$zs,
           main = paste0(Y_name, " module ", c, ", acceptance curve"),
           xlab = "Index", ylab = "z-statistic")
      points(results$commzs$Y_zs[[c]]$zs[1:YL], col = "red")
      dev.off()

    }

  }
  
  # Doing best match
  
  dX <- ncol(X)
  dY <- ncol(Y)
  
  # Preparing component list for scoring
  bg_list <- setdiff(1:(dX + dY), unlist(component_list))
                          
  
  # Preparing BMD results for scoring
  results_bg <- c(results$background$X_bg, results$background$Y_bg)
  results_comms <- lapply(seq_along(results$communities[[1]]), 
                             function (c) {
                               c(results$communities$X_sets[[c]],
                                 results$communities$Y_sets[[c]])
                             }
  )
  
  BMD_bestmatch <- best_match_bimodule(component_list, results_comms, 
                                       bg_list, results_bg)
  
  return(BMD_bestmatch)

}






