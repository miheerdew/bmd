library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(bmdCpp)

#sourceCpp("correlation.cpp")

cbceNW <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, Cpp = TRUE, verbose = TRUE, generalOutput = TRUE,
                    updateOutput = TRUE, throwInitial = TRUE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                    updateMethod = 1, inv.length = TRUE, add_rate = 1, start_nodes = NULL,
                    calc_full_cor=FALSE, loop_limit = Inf, parallel = FALSE, twoSided = FALSE) {
  
  if (FALSE) {
    alpha = 0.05
    OL_thres = 0.9
    tag = NULL
    verbose = TRUE
    generalOutput = TRUE
    OL_tol = Inf
    Cpp = TRUE
    Dud_tol = Inf
    time_limit = Inf
    updateMethod = 1
    updateOutput = TRUE
    throwInitial = TRUE
    inv.length = TRUE
    time_limit = Inf
    loop_limit = Inf
    add_rate = 1
    calc_full_cor = TRUE
    parallel = FALSE
    start_nodes = NULL
    twoSided = FALSE
  }
  
  #-----------------------------------------------------------------------------
  # Setup 
  
  start_second <- proc.time()[3]
  td <- tempdir()
  cat("#-------------------\n")
  
  if (generalOutput)
    cat("Setup\n")
  
  X <- scale(X); Y <- scale(Y)
  X3 <- X^3; X2 <- X^2; X4ColSum <- colSums(X^4)
  Y3 <- Y^3; Y2 <- Y^2; Y4ColSum <- colSums(Y^4)
  
  dx <- ncol(X)
  dy <- ncol(Y)
  n  <- nrow(X)
  
  Xindx <- 1:dx
  Yindx <- (dx + 1):(dx + dy)
  
  if(calc_full_cor){
    if (generalOutput)
      cat("Calculating full cross correlation matrix\n")
    full_xy_cor <- cor(X, Y)
  }
  
  source("auxiliary.R", local = TRUE)
  source("bh_reject.R", local = TRUE)
  source("pvals.R", local = TRUE)
  source("initialize.R", local = TRUE)
  source("extractNW.R", local = TRUE)
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  if (generalOutput)
    cat("Beginning method.\n\n")
  
  # Getting node orders.
  Ysum <- Y %*% rep(1,dy) / dy
  Xsum <- X %*% rep(1,dx) / dx
  cor_X_to_Ysums <- as.vector(t(Ysum) %*% X)
  cor_Y_to_Xsums <- as.vector(t(Xsum) %*% Y)
  
  if (!twoSided) {
    extractord <- c(Xindx, Yindx)[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                          decreasing = TRUE)]
  } else {
    extractord <- sample(c(Xindx, Yindx))
  }
  
  if (!is.null(start_nodes))
    extractord <- extractord[extractord %in% start_nodes]
  
  # Initializing control variables
  stop_extracting <- FALSE
  OL_count <- 0
  Dud_count <- 0
  clustered <- integer(0)
  comms <- NULL
  
  # Extracting
  if (parallel) {
    ticp <- proc.time()[3]
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, library(bmdCpp))
    registerDoParallel(cl)
    extract_res <- foreach(i = extractord) %dopar% {
      extract(i, print_output = FALSE, interact = TRUE)
    }
    stopCluster(cl)
    tocp <- proc.time()[3]
  } else {
    tic <- proc.time()[3]
    extract_res <- lapply(extractord, extract, interact = TRUE)
    extract_res <- extract_res[order(extractord)]
    toc <- proc.time()[3]
  }
  
  #-----------------------------------------------------------------------------
  # Clean-up and return -------------------------------------------------------
  
  if (generalOutput)
    cat("Cleaning up.\n")
  
  # Getting final sets and making report
  final.sets <- lapply(extract_res, function (R) R$StableComm)
  endtime <- proc.time()[3]
  report <- list(OL_count, Dud_count, endtime - start_second)
  names(report) <- c("OL_count", "Dud_count", "timer")
  
  # Removing blanks and trivial sets
  nonNullIndxs <- which(unlist(lapply(final.sets, length)) > 0)
  if (length(nonNullIndxs) == 0) {
    returnList <- list("communities" = list("X_sets" = NULL,
                                            "Y_sets" = NULL),
                       "background" = 1:(dx + dy),
                       "extract_res" = extract_res,
                       "finalIndxs" = integer(0),
                       "final.sets" = final.sets,
                       "report" = report)
    return(returnList)
  }
  nonNullComms <- final.sets[nonNullIndxs]
  OLfilt <- filter_overlap(nonNullComms, tau = OL_thres, inv.length = inv.length)
  finalComms <- OLfilt$final_comms
  finalIndxs <- nonNullIndxs[OLfilt$kept_comms]
  
  returnList <- list("extract_res" = extract_res,
                     "finalIndxs" = finalIndxs,
                     "report" = report)
  
  cat("#-------------------\n\n\n")
  
  return(returnList)
  
}
