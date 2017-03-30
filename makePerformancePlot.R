tonmi <- function(x, a) {
  (1 / (1 - x + a) - 1 / (1 + a)) / (1 / a - 1 / (1 + a)) 
}

makePerformancePlot <-  function (fn = NULL, meanMat, sdMat = NULL,
                                  legPos = "topright", xvals, tnmi = FALSE,
                                  tnmi_a = 0.05, pns = plot_names,
                                  legCex = 2, pchs = 1:nrow(meanMat),
                                  lwd = 1, legLwd = 2, yRange = NULL,
                                  cex = 1, new = TRUE, doLegend = TRUE,
                                  xRange, plotFile = TRUE,
                                  ... ) {
  
  # Calculating plot limits
  
  if (!new) {
    upperLimit <- max(meanMat + sdMat)
    lowerLimit <- min(meanMat - sdMat)
  } else {
    upperLimit <- max(meanMat)
    lowerLimit <- min(meanMat)
  }
  if (plotFile) {
    png(fn, width = 1000, height = 1000)
    
    par(oma= rep(0,4),
        mar=c(12,12,6,4),
        mgp=c(8,3,0),
        mfrow = c(1,1))
  }
  
  if (is.null(yRange)) {
    yRange <- c(lowerLimit, upperLimit)
  }
  
  if (tnmi) {
    yRange <- c(0, 1)
    meanMat <- tonmi(meanMat, tnmi_a)
  }
  
  plot(0.5, 0.5, col="white", pch='.',
       xlim = xRange, ylim = yRange, ...)
  
  for(j in 1:length(pns)){
    
    if (new) {
      
      meth = pns[j]
      points(xvals, meanMat[meth, ], pch = pchs[j], cex = cex, col = colPal[j])
      lines(xvals, meanMat[meth, ], lty = j + 1, lwd = lwd, col = colPal[j])
      if (tnmi) {
        abline(h = tonmi(0.95, tnmi_a), lwd = lwd)
        text(y = tonmi(0.95, tnmi_a), pos = 1, cex = cex.lab,
             x = paramVec[round(length(paramVec) * 0.75)],
             labels = c("0.95"))
      }
      
    } else {
    
      meth = pns[j]
      plotCI(x = xvals,
             y = meanMat[meth, ],
             uiw = sdMat[meth, ],         
             pch = pchs[j],
             pt.bg="black",
             cex=cex,
             lty=1,
             lwd = lwd,
             gap=0,
             type="o",
             sfrac=0.005,
             add=TRUE,
             col = colPal[j])
      
    }
    
  }
  
  if (doLegend) {
    if (new) {
      
      legend(x = legPos, legend = pns, lty = 1 + 1:length(pns),
             cex = legCex, lwd = legLwd, pch = pchs, col = colPal)
      
    } else {
    
    legend(x = legPos,
           legend = pns,
           col = colPal,
           lty = 1,
           lwd = lwd,
           pt.cex = cex,
           pch = 1:length(pns),
           cex = legCex)
      
    }
  }
  
  if (plotFile)
    dev.off()
  
  return(NULL)
  
}