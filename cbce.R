library(Rcpp)
library(RcppParallel)
library(foreach)
library(doParallel)
library(bmdCpp)
source("auxiliary.R")
source("bh_reject.R")
source("pvals.R")
source("initialize.R")
source("extract.R")

#sourceCpp("correlation.cpp")

bmd <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, Cpp = TRUE, verbose = TRUE, generalOutput = TRUE,
                 updateOutput = TRUE, throwInitial = TRUE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                 updateMethod = 1, inv.length = TRUE, add_rate = 1, start_nodes = NULL,
                 calc_full_cor=FALSE, loop_limit = Inf, parallel = FALSE) {
  
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
  
  #-------------------------------------------------------------------------------
  # Extractions
  
  if (generalOutput)
    cat("Beginning method.\n\n")
  
  # Getting node orders.
  Ysum <- Y_scaled %*% rep(1,dy) / dy
  Xsum <- X_scaled %*% rep(1,dx) / dx
  cor_X_to_Ysums <- as.vector(t(Ysum) %*% X_scaled)
  cor_Y_to_Xsums <- as.vector(t(Xsum) %*% Y_scaled)
  
  extractord <- c(Xindx, Yindx)[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                      decreasing = TRUE)]
  
  if (!is.null(start_nodes))
    extractord <- extractord[extractord %in% start_nodes]
  
  # Initializing control variables
  stop_fn <- file.path(td, "stop_extracting.txt")
  OL_fn <- file.path(td, "OL_count.txt")
  Dud_fn <- file.path(td, "Dud_count.txt")
  comm_dn <- file.path(td, "comm_dn")
  file.create(stop_fn, OL_fn, Dud_fn)
  dir.create(comm_dn, showWarnings = FALSE)
  writeLines("0", OL_fn); writeLines("0", Dud_fn); writeLines("FALSE", stop_fn)
  
  # Extracting
  if (parallel) {
    ticp <- proc.time()[3]
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
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
  
  # Getting final sets and counts
  final.sets <- lapply(extract_res, function (R) R$StableComm)
  OL_count <- as.numeric(readLines(OL_fn))
  Dud_count <- as.numeric(readLines(Dud_fn))
  
  # Removing temp files
  file.remove(stop_fn, OL_fn, Dud_fn)
  unlink(comm_dn, recursive = TRUE)
  
  # Making report
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
