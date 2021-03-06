library(doParallel)
library(foreach)
source("makeVars.R")
source("stdize.R")
Rcpp::sourceCpp("bmd_helper.cpp")

bmd_cpp <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, cp_cor = TRUE, verbose = TRUE,
                     updateOutput = TRUE, throwInitial = TRUE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                     updateMethod = 5, initializeMethod = 3, inv.length = TRUE, add_rate = 1,
                     calc_full_cor=FALSE, loop_limit = Inf, parallel = FALSE, conserv=TRUE) {

  if (FALSE) {
    alpha = 0.05
    OL_thres = 0.9
    tag = NULL
    verbose = TRUE
    OL_tol = Inf
    Dud_tol = Inf
    time_limit = Inf
    updateMethod = 5
    initializeMethod = 3
    updateOutput = TRUE
    throwInitial = TRUE
    inv.length = TRUE
    time_limit = Inf
    loop_limit = Inf
    add_rate = 1
    calc_full_cor=TRUE
    parallel = FALSE
    conserv = TRUE
  }

  start_second <- proc.time()[3]


  #-------------------------------------------------------------------------------
  # Auxiliary Functions

  symdiff <- function (s1, s2) {
    return(union(setdiff(s1, s2), setdiff(s2, s1)))
  }

  jaccard <- function (s1, s2) {
    return(length(symdiff(s1, s2)) / length(union(s1, s2)))
  }

  filter_overlap <- function (comms, tau, inv.length = FALSE) {

    K <- length(comms)
    if (inv.length) {
      scores <- 1 / unlist(lapply(comms, length))
    } else {
      scores <- unlist(lapply(comms, length))
    }


    jaccard_mat0 <- matrix(0, K, K)
    for (i in 1:K) {
      for (j in 1:K) {
        jaccard_mat0[i, j] <- length(intersect(comms[[i]], comms[[j]])) /
          length(comms[[i]])
      }
    }

    jaccard_mat <- jaccard_mat0
    diag(jaccard_mat) <- 0
    max_jacc <- max(jaccard_mat)
    deleted_comms <- integer(0)

    while (max_jacc > tau) {

      inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]

      # keep comm with larger score
      delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
      jaccard_mat[delete_comm, ] <- 0
      jaccard_mat[, delete_comm] <- 0
      deleted_comms <- c(deleted_comms, delete_comm)
      max_jacc <- max(jaccard_mat)

    }

    kept_comms <- setdiff(1:K, deleted_comms)

    return(list("final_comms" = comms[kept_comms],
                "kept_comms" = kept_comms))

  }


  #-----------------------------------------------------------------------------
  # Setup calculations & update functions

  X_orig <- X; Y_orig <- Y;

  cat("Setting up calculations\n")

  X <- scale(X); Y <- scale(Y)
  X3 <- X^3; X2 <- X^2; X4ColSum <- colSums(X^4)
  Y3 <- Y^3; Y2 <- Y^2; Y4ColSum <- colSums(Y^4)

  if(calc_full_cor){
    cat("Calculating full cross correlation matrix\n")
    full_xy_cor <- cor(X, Y)
  }

  dx <- ncol(X)
  dy <- ncol(Y)
  n  <- nrow(X)

  Xindx <- 1:dx
  Yindx <- (dx + 1):(dx + dy)

  initialize <- function(u){
    if(u > dx){
      #Test X
      u <- u - dx
      return(initializeC(n,
                        if(calc_full_cor) full_xy_cor[,u] else cor(X, Y[,u]),
                        alpha,
                        conserv))
    } else {
      #Test Y
      return(initializeC(n,
                        if(calc_full_cor) t(full_xy_cor[u,]) else cor(Y, X[,u]),
                        alpha,
                        conserv) + dx)
    }
  }

  if(calc_full_cor){
    update5 <- function(A, B=NULL){
      if(min(A) > dx){
        #Test X
        A <- A - dx
        return(updateC(X, Y[,A,drop=FALSE], X4ColSum, X2, X3,
                       full_xy_cor[,A, drop=FALSE],
                       alpha, conserv))
      } else {
        #Test Y
        return(updateC(Y, X[,A,drop=FALSE], Y4ColSum, Y2, Y3,
                       t(full_xy_cor[A,,drop=FALSE]),
                       alpha, conserv) + dx)
      }
    }
  } else {
    update5 <- function(A, B=NULL){
      if(min(A) > dx){
        #Test X
        A <- A - dx
        return(updateC(X, Y[,A,drop=FALSE], X4ColSum, X2, X3,
                       cor(X, Y[,A,drop=FALSE]),
                       alpha, conserv))
      } else {
        #Test Y
        return(updateC(Y, X[,A,drop=FALSE], Y4ColSum, Y2, Y3,
                       cor(Y, X[,A,drop=FALSE]),
                       alpha, conserv) + dx)
      }
    }
  }

  pvals <- function(A){
    if(min(A) > dx){
      #Test X
      A <- A - dx
      return(pvalsC(X, Y[,A,drop=FALSE], X4ColSum, X2, X3,
              if(calc_full_cor) full_xy_cor[,A,drop=FALSE] else cor(X, Y[,A])))
    } else {
      #Test Y
      return(pvalsC(Y, X[,A,drop=FALSE], Y4ColSum, Y2, Y3,
              if(calc_full_cor) t(full_xy_cor[A,drop=FALSE]) else cor(Y, X[,A])))
    }
  }

  #-----------------------------------------------------------------------------
  # Extract function

  extract <- function (indx, interact = FALSE, print_output = verbose) {

    if (interact && (indx %in% clustered || stop_extracting)) return(integer(0))

    if (print_output && interact) {
      cat("\n#-----------------------------------------\n\n")
      cat("extraction", which(extractord == indx), "of", length(extractord), "\n\n")
      cat("OL_count = ", OL_count, "\n")
      cat("Dud_count = ", Dud_count, "\n")
      cat(paste0(sum(clustered <= dx), " X vertices in communities.\n"))
      cat(paste0(sum(clustered > dx), " Y vertices in communities.\n"))
    }

    B0x <- initialize(indx)
    if (length(B0x) <= 1) {
      if (interact) Dud_count <<- Dud_count + 1
      return(integer(0))
    }
    B0y <- update5(B0x, indx)

    # Initializing extraction loop
    B_oldx <- B0x; B_oldy <- B0y
    B_old <- c(B_oldx, B_oldy)
    initial_set <- B_old
    B_new <- c(Xindx, Yindx)
    chain <- list(B_old)
    consec_jaccards <- mean_jaccards <- NULL
    consec_sizes <- list(c(length(B_oldx), length(B_oldy)))
    found_cycle <- found_break <- cycledSets <- NULL
    did_it_cycle <- FALSE
    itCount <- 0

    repeat {

      itCount <- itCount + 1

      if (updateMethod != 7) {
        B_newx <- update5(B_oldy, B_oldx)
        if (length(B_newx) >= 1) {
          B_newy <- update5(B_newx, B_oldy)
        } else {
          B_newy <- integer(0)
          break
        }
        B_new <- c(B_newx, B_newy)

        if (length(B_newy) == 0)
          break
      } else {

        if (comm_indx > dx) {
          Xpvals <- pvals(B_oldy)
          Ypvals <- pvals(B_oldx)
          B_new <- bh_rejectC(c(Xpvals, Ypvals), alpha, conserv)
          B_newx <- B_new[B_new <= dx]
          B_newy <- B_new[B_new > dx]
        } else {
          Xpvals <- pvals(B_oldx)
          Ypvals <- pvals(B_oldy)
          B_new <- bh_rejectC(c(Xpvals, Ypvals), alpha, conserv)
          B_newx <- B_new[B_new > dx]
          B_newy <- B_new[B_new <= dx]
        }

        if (length(B_newy) * length(B_newx) == 0) {
          B_newy <- B_newx <- integer(0)
          break
        }

      }

      consec_jaccard <- jaccard(B_new, B_old)
      jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
      found_cycle <- c(found_cycle, FALSE)
      found_break <- c(found_break, FALSE)
      consec_jaccards <- c(consec_jaccards, consec_jaccard)
      mean_jaccards <- c(mean_jaccards, mean(jaccards))
      consec_sizes <- c(consec_sizes, list(c(length(B_newx), length(B_newy))))

      # Checking for cycles (4.4.1 in CCME paper)
      if (jaccard(B_new, B_old) > 0) { # Otherwise loop will end naturally

        if (sum(jaccards == 0) > 0) { # Cycle has been found

          found_cycle[itCount] <- TRUE
          did_it_cycle <- TRUE

          if (updateOutput)
            cat("---- Cycle found")

          Start <- max(which(jaccards == 0))
          cycle_chain <- chain[Start:length(chain)]

          # Checking for cycle break (4.4.1a)
          seq_pair_jaccards <- rep(0, length(cycle_chain) - 1)
          for (j in seq_along(seq_pair_jaccards)) {
            seq_pair_jaccards[j] <- jaccard(chain[[Start + j - 1]], chain[[Start + j]])
          }
          if (sum(seq_pair_jaccards > 0.5) > 0) {# then break needed
            cat(" ---- Break found\n")
            found_break[itCount] <- TRUE
            B_new <- integer(0)
            break
          }

          # Create conglomerate set (and check, 4.4.1b)
          B_J <- unique(unlist(cycle_chain))
          B_J_check <- unlist(lapply(chain, function (B) jaccard(B_J, B)))
          if (sum(B_J_check == 0) > 0) {
            cat(" ---- Old cycle\n")
            break
          } else {
            cat(" ---- New cycle\n")
            B_oldx <- B_J[B_J <= dx]
            B_oldy <- B_J[B_J > dx]
          }

        } else {
          # From checking jaccards to cycle_chain; if not, then can set B_oldx
          # and B_oldy to the update and restart.
          B_oldx <- B_newx
          B_oldy <- B_newy
        }

      } else { # From checking B_new to B_old; if not, B_new = B_old and:
        break
      }

      B_old <- c(B_oldx, B_oldy)
      chain <- c(chain, list(B_new))

    }

    # Storing B_new and collecting update info
    if (length(B_newx) * length(B_newy) == 0)
      B_new <- integer(0)
    if (interact) final.sets[[indx]] <<- B_new
    update_info <- list("mean_jaccards" = mean_jaccards,
                        "consec_jaccards" = consec_jaccards,
                        "consec_sizes"= consec_sizes,
                        "found_cycle" = found_cycle,
                        "found_break" = found_break)
    if (interact) clustered <<- union(clustered, B_new)

    # Checking overlap with previous sets
    if (interact && sum(plugged) > 0) {
      OL_check <- unlist(lapply(final.sets[plugged],
                                function (C) jaccard(B_new, C)))
      if (sum(OL_check < 1 - OL_thres) > 0)
        OL_count <<- OL_count + 1
    }

    # Noting which final.sets are filled; doing checks
    if (interact) plugged[indx] <<- length(B_new) > 0
    current_time <- proc.time()[3] - start_second
    if (interact && current_time > time_limit || Dud_count > Dud_tol || OL_count > OL_tol)
      stop_extracting <<- TRUE

    return(list("StableComm" = B_new,
                "update_info" = update_info,
                "initial_set" = initial_set,
                "itCount" = itCount, "did_it_cycle" = did_it_cycle,
                "current_time" = current_time))
  }

  #-------------------------------------------------------------------------------
  # Extractions

  cat("Beginning method.\n\n")

  # Getting node orders.
  Ysum <- Y %*% rep(1,dy) / dy
  Xsum <- X %*% rep(1,dx) / dx
  cor_X_to_Ysums <- as.vector(t(Ysum) %*% X)
  cor_Y_to_Xsums <- as.vector(t(Xsum) %*% Y)

  extractord <- c(Xindx, Yindx)[order(c(cor_X_to_Ysums, cor_Y_to_Xsums),
                                      decreasing = TRUE)]

  # Initializing control variables
  clustered <- integer(0)
  stop_extracting <- FALSE
  plugged <- logical(dx + dy)
  OL_count <- 0
  Dud_count <- 0
  final.sets <- rep(list(integer(0)), length(extractord))

  # Extracting
  if (parallel) {
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    extract_res <- foreach(i = extractord, .combine='c') %dopar% {
      extract(i, interact = TRUE, print_output = FALSE)
    }
    stopCluster(cl)
  } else {
    extract_res <- lapply(extractord, extract, interact = TRUE)
    extract_res <- extract_res[order(extractord)]
  }

  #-----------------------------------------------------------------------------
  #  Clean-up and return -------------------------------------------------------

  endtime <- proc.time()[3]
  report <- list(OL_count, Dud_count, endtime - start_second)
  names(report) <- c("OL_count", "Dud_count", "timer")

  # Removing blanks and trivial sets
  cat("removing overlap and unpacking comms...\n")
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

  X_sets <- lapply(finalComms, function (C) C[C <= dx])
  Y_sets <- lapply(finalComms, function (C) C[C > dx])
  X_bg <- setdiff(Xindx, unlist(X_sets))
  Y_bg <- setdiff(Yindx, unlist(Y_sets))

  returnList <- list("communities" = list("X_sets" = X_sets,
                                          "Y_sets" = Y_sets),
                     "background" = list("X_bg" = X_bg,
                                         "Y_bg" = Y_bg),
                     "extract_res" = extract_res[finalIndxs],
                     "finalIndxs" = finalIndxs,
                     "final.sets" = final.sets,
                     "report" = report)
  return(returnList)

}
