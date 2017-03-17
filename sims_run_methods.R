library(MASS)
library(Matrix)
library(lpbrim)
source("bmd.R")
source("formatBRIM.R")
source("sims_config.R")
source("sims_debug_tools.R")

for (n in ns) {

  # Loading data
  load(dataset_fname(n))

  if (doBRIM) {

    # Running lpbrim
    BRIMtime <- proc.time()[3]

    # Making large correlation matrix
    n <- nrow(X)
    p1 <- ncol(X)
    p2 <- ncol(Y)
    p <- p1 + p2
    crossCors <- cor(X, Y)
    crossTstats <- crossCors * sqrt(n - 2) / sqrt(1 - crossCors^2)
    crossPvals <- as.vector(pt(crossTstats, n - 2, lower.tail = FALSE))
    crossPvals_adj <- crossPvals * p1 * p2 / rank(crossPvals)
    crossSigs <- as.integer(crossPvals_adj <= 0.1)
    crossSigMat <- matrix(crossSigs, ncol = p2)
    largeM <- matrix(integer(p^2), ncol = p)
    largeM[1:p1, (p1 + 1):p] <- crossSigMat
    largeM[(p1 + 1):p, 1:p1] <- t(crossSigMat)

    # Making sure no rowSums are zero
    zeroDegs <- rowSums(largeM) == 0
    largeM_red <- largeM[!zeroDegs, !zeroDegs]
    BRIMmod <- findModules(largeM_red, iter = 5, sparse = FALSE)
    BRIMformat <- formatBRIM(BRIMmod, zeroDegs, ncol(X))
    BRIMresults1 <- list("communities" = BRIMformat[[1]])
    BRIMresults2 <- list("communities" = BRIMformat[[2]])
    BRIMtime <- proc.time()[3] - BRIMtime

    save(BRIMtime, BRIMresults1, BRIMresults2,
         file = results_fname(n, method="brim"))

  }

  # Running bmd
  BMDtime <- proc.time()[3]
  suppressWarnings(
  BMDresults <- bmd(X, Y, tag = n,
                    saveDir = file.path(saveDir, "BMD_saves"),
                    updateMethod = 5, initializeMethod = 3,
                    Dud_tol = 10, OL_tol = 10, time_limit = 1800,
                    bmd_index = bmd_index, calc_full_cor=TRUE)
  )

  BMDtime <- proc.time()[3] - BMDtime
  save(BMDtime, BMDresults,
       file = results_fname(n, method="bmd"))
  
  if (doKB) {
    BMDtime <- proc.time()[3]
    suppressWarnings(
    BMDresults <- bmd(X, Y, tag = n,
                      saveDir = file.path(saveDir, "BMD_saves"),
                      updateMethod = 6, initializeMethod = 3,
                      Dud_tol = 10, OL_tol = 10, time_limit = 1800,
                      bmd_index = bmd_index, calc_full_cor=TRUE)
    )
    BMDtime <- proc.time()[3] - BMDtime
    save(BMDtime, BMDresults,
         file = results_fname(n, method="bmd-kb"))
  }
  
}
