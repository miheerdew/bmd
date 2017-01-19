library(gplots)
library(RColorBrewer)
### Brewing colors for plotting
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

plotLargeMat <- function (Mat, X_sets, Y_sets, maxPlot = 1e4,
                           col = r, sepcol = "black", ...) {
  
  
  X_grabs <- unlist(X_sets)
  Y_grabs <- unlist(Y_sets)
  nX <- length(X_grabs)
  nY <- length(Y_grabs)
  
  Xsizes <- unlist(lapply(X_sets, length))
  Ysizes <- unlist(lapply(Y_sets, length))

  # Determining shrinkage factors
  if (nY * nX > 10 * maxPlot) {
    
    c <- sqrt(maxPlot / (nY * nX))
    chunkSize <- ceiling(1 / c)
    
    # Dealing with X
    nXchunks <- floor(nX / chunkSize)
    nXremove <- nX - nXchunks * chunkSize
    if (nXremove > 0) {
      largest <- which.max(Xsizes)
      to_remove <- sample(Xsizes[largest], nXremove, replace = FALSE)
      X_sets[[largest]] <- X_sets[[largest]][-to_remove]
    }
    
    # Dealing with Y
    nYchunks <- floor(nY / chunkSize)
    nYremove <- nY - nYchunks * chunkSize
    if (nYremove > 0) {
      largest <- which.max(Ysizes)
      to_remove <- sample(Ysizes[largest], nYremove, replace = FALSE)
      Y_sets[[largest]] <- Y_sets[[largest]][-to_remove]
    }
    
    # Re-doing counting
    X_grabs <- unlist(X_sets)
    Y_grabs <- unlist(Y_sets)
    nX <- length(X_grabs)
    nY <- length(Y_grabs)
    
    Xsizes <- unlist(lapply(X_sets, length))
    Ysizes <- unlist(lapply(Y_sets, length))
    
    Xsep <- cumsum(unlist(lapply(X_sets, length)))
    Ysep <- cumsum(unlist(lapply(Y_sets, length)))
    Xsep <- Xsep[-length(Xsep)]
    Ysep <- Ysep[-length(Ysep)]
    
    # Making 2nd crossCors
    cat("Making largeMat...\n")
    largeMat <- Mat[X_grabs, Y_grabs]
    largeMat_red1 <- matrix(0, nXchunks, nY)
    
    # Filling largeMat_red1
    cat("Filling largeMat_red1...\n")
    pos <- 1
    chunkCount <- 1
    while (pos < nX) {
      indxs <- pos:(pos + chunkSize - 1)
      largeMat_red1[chunkCount, ] <- colSums(largeMat[indxs, ])
      pos <- pos + chunkSize
      chunkCount <- chunkCount + 1
    }
    
    # Filling largeMat_red2
    largeMat_red2 <- matrix(0, nXchunks, nYchunks)
    pos <- 1
    chunkCount <- 1
    while (pos < nY) {
      indxs <- pos:(pos + chunkSize - 1)
      largeMat_red2[ , chunkCount] <- rowSums(largeMat_red1[ , indxs])
      pos <- pos + chunkSize
      chunkCount <- chunkCount + 1
    }
    
    # Averaging
    largeMat_red2 <- largeMat_red2 / chunkSize^2
    
    # Re-doing sep points
    nX_red <- nrow(largeMat_red2)
    nY_red <- ncol(largeMat_red2)
    cX <- nX_red / nX
    cY <- nY_red / nY
    Xsep_red <- round(Xsep * cX)
    Ysep_red <- round(Ysep * cY)
    Xsep_red <- pmin(Xsep_red, nX_red)
    Xsep_red <- pmax(Xsep_red, 0)
    Ysep_red <- pmin(Ysep_red, nY_red)
    Ysep_red <- pmax(Ysep_red, 0)
    
    Xsep <- Xsep_red
    Ysep <- Ysep_red
    rm(largeMat, largeMat_red1)
    gc()
    
  } else {
    largeMat_red2 <- Mat[X_grabs, Y_grabs]
    Xsep <- cumsum(unlist(lapply(X_sets, length)))
    Ysep <- cumsum(unlist(lapply(Y_sets, length)))
    Xsep_red <- Xsep[-length(Xsep)]
    Ysep_red <- Ysep[-length(Ysep)]
  }
  
  heatmap.2(largeMat_red2, Rowv = FALSE, Colv = FALSE, dendrogram = "none",
            trace = "none", labRow = "", labCol = "",
            colsep = Ysep_red, rowsep = Xsep_red, col = col, sepcol = sepcol, ...)
  
}