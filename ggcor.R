library(reshape2)
library(ggplot2)

ggcor <- function (cormat, fn, removeLowerTri = FALSE, 
                   fisher = TRUE, n = NULL, title = "corrstats") {
  
  diag(cormat) <- 0
  
  if (removeLowerTri)
    cormat[lower.tri(cormat)] <- NA
  
  # Melt the correlation matrix
  melted_cormat <- melt(cormat, na.rm = TRUE)
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
  
  
  
  png(fn)
  print(ggplot(data = melted_cormat, aes(x = Var2, y = Var1, fill = value)) + 
        geom_tile() + colorscale + ggtitle(title)
        
  )
  dev.off()
  
}