library(reshape2)
library(ggplot2)

ggcor <- function (cormat, fn = NULL, removeLowerTri = FALSE, 
                   fisher = TRUE, n = NULL, title = "corrstats") {
  
  diag(cormat) <- 0
  
  # Shrinking the matrix if necessary
  if (nrow(cormat) == ncol(cormat) && nrow(cormat) > 1000) {
    cat("filling smaller cormat...\n")
    blocksize <- ceiling(nrow(cormat) / 1000)
    newcormat <- matrix(0, 1000, 1000)
    mposi <- 1
    for (i in 1:1000) {
      if (mposi > nrow(cormat)) break
      if (i %% 100 == 0)
        cat("i =", i, "/ 1000\n")
      mposj <- 1
      for (j in 1:1000) {
        if (mposj > nrow(cormat)) break
        indxi <- mposi:min(nrow(cormat), (mposi + blocksize - 1))
        indxj <- mposj:min(nrow(cormat), (mposj + blocksize - 1))
        newcormat[i, j] <- mean(cormat[indxi, indxj])
        mposj <- mposj + blocksize
      }
      mposi <- mposi + blocksize
    }
  }
  
  if (removeLowerTri)
    newcormat[lower.tri(cormat)] <- NA
  
  # Melt the correlation matrix
  melted_cormat <- melt(newcormat, na.rm = TRUE)
  melted_cormat$Var1 <- max(melted_cormat$Var1) - melted_cormat$Var1 + 1
  
  if (fisher) {
    if (is.null(n))
      stop('n must be supplied to calculate scaled fisher values\n')
    melted_cormat$value <- atanh(melted_cormat$value) * sqrt(n - 3)
    maxabs <- ceiling(max(abs(melted_cormat$value)))
    plotrange <- c(-maxabs, maxabs)
    colorscale <- scale_fill_gradient2(low = "blue", 
                                       high = "red", 
                                       mid = "white", 
                                       midpoint = 0, limits = plotrange,
                                       space = "Lab", 
                                       name = "Fisher")
  } else {
    
    maxabs <- max(abs(melted_cormat$value))
    plotrange <- c(-maxabs, maxabs)
    colorscale <- scale_fill_gradient2(low = "blue", 
                                       high = "red", 
                                       mid = "white", 
                                       midpoint = 0, limits = plotrange,
                                       space = "Lab", 
                                       name = "Corr")
    
  }
  
  p <- ggplot(data = melted_cormat, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + colorscale + ggtitle(title)
  if (!is.null(fn)) {
    png(fn)
    print(p)
    dev.off()
  } else {
    return(p)
  }
}