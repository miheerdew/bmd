library(Rcpp)
library(RcppParallel)
source("makeVars.R")
source("stdize.R")
#sourceCpp("correlation.cpp")

bmd <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, saveDir = NULL, cp_cor = TRUE,
                 updateOutput = TRUE, throwInitial = TRUE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                 updateMethod = 5, initializeMethod = 3, inv.length = TRUE, add_rate = 1,
                 bmd_index=NULL, calc_full_cor=FALSE, loop_limit = Inf) {
  # bmd_index : A function that maps each vertex to the index of the bimodule
  #           it contains.
  
  if (FALSE) {
    alpha = 0.05
    OL_thres = 0.9
    tag = NULL
    saveDir = NULL
    OL_tol = Inf
    Dud_tol = Inf
    time_limit = 18000
    updateMethod = 5
    initializeMethod = 3
    updateOutput = TRUE
    throwInitial = TRUE
    inv.length = TRUE
    time_limit = Inf
    loop_limit = Inf
    add_rate = 1
    calc_full_cor=TRUE
    bmd_index=NULL
  }
  
  start_second <- proc.time()[3]
  starttime <- Sys.time()
  if (!is.null(saveDir) && !dir.exists(saveDir))
    dir.create(saveDir)
  
  
  #-------------------------------------------------------------------------------
  # Auxiliary Functions --------------------------------------------------------
  
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
  
  bh_reject <- function (pvals, alpha, conserv = TRUE) {
    
    m <- length(pvals)
    
    if (!conserv) {
      pvals_adj <- m * pvals / rank(pvals)
    } else {
      mults <- sum(1 / c(1:m))
      pvals_adj <- mults * m * pvals / rank(pvals)
    }
    
    if (sum(pvals_adj <= alpha) > 0) {
      thres <- max(pvals[pvals_adj <= alpha])
      return(which(pvals <= thres))
    } else {
      return(integer(0))
    }
  }
  
  # Setup calculations    
  
  n <- nrow(X)
  dx <- ncol(X)
  dy <- ncol(Y)
  X_scaled <- scale(X)
  Y_scaled <- scale(Y)
  
  if(calc_full_cor || dx * dy <= 1e7){
    cat("Calculating the full cross correlation matrix.\n")
    full_xy_cor = cor(X,Y)
  }
  
  Xindx <- 1:dx
  Yindx <- (dx + 1):(dx + dy)
  
  cross_cor <- function(X_indx=1:dx, Y_indx=1:dy){
    if (calc_full_cor){
      return(full_xy_cor[X_indx, Y_indx, drop = FALSE])
    } else {
      return(cor(as.matrix(X[,X_indx, drop = FALSE]), as.matrix(Y[,Y_indx, drop = FALSE])))
    }
    
  }
  
  update5 <- function (B, A = NULL, justpvals = FALSE) {
    
    if (length(B) == 0)
      return(integer(0))
    
    test_X <- min(B) > dx
    nFixd <- length(B)
    
    if (test_X) {
      
      # Getting fixed matrix
      fixdIndx <- match(B, Yindx)
      fixdMat <- Y_scaled[ , fixdIndx, drop = FALSE]
      
      # Calculating the variances
      {
        # General calcs
        xyCors <- cross_cor(Y_indx = fixdIndx)
        #if (calc_full_cor){
        #  xyCors <- full_xy_cor[ , fixdIndx, drop = FALSE]
        #} else {
        #  xyCors <- cor(X, Y[,fixedIndx, drop = FALSE])
        #}
        y4 <- colSums(X_scaled^4)
        xRowSum <- rowSums(fixdMat)
        xRowSum2 <- tcrossprod(xyCors, fixdMat^2)
        
        # Calc for star 1
        star1 <- crossprod(X_scaled^2, xRowSum^2)
        
        # Calc for star 2
        star2 <- y4 * rowSums(xyCors)^2
        
        # Calc for star 3
        star3 <- 2 * rowSums(xyCors) * colSums(X_scaled^2 * t(xRowSum2))
        
        # Calc for star 4
        star4 <- rowSums(xRowSum2^2)
        
        # Calc for dagger 1
        dagger1 <- rowSums(xyCors) * crossprod(X_scaled^3, xRowSum)
        
        # Calc for dagger 2
        dagger2 <- colSums(xRowSum * t(xRowSum2) * X_scaled)
      }
      
      
    } else {
      
      # Getting indices
      fixdIndx <- match(B, Xindx)
      fixdMat <- X_scaled[ , fixdIndx, drop = FALSE]
      
      # Calculating the variances
      {
        # General calcs
        xyCors <- t(cross_cor(X_indx=fixdIndx))
        #if (calc_full_cor){
        #  xyCors <- full_xy_cor[fixdIndx, , drop = FALSE]
        #} else {
        #  xyCors <- cor(X[,fixdIndx, drop = FALSE], Y)
        #}
        #xyCors <- t(xyCors)
        y4 <- colSums(Y_scaled^4)
        xRowSum <- rowSums(fixdMat)
        xRowSum2 <- tcrossprod(xyCors, fixdMat^2)
        
        # Calc for star 1
        star1 <- crossprod(Y_scaled^2, xRowSum^2)
        
        # Calc for star 2
        star2 <- y4 * rowSums(xyCors)^2
        
        # Calc for star 3
        star3 <- 2 * rowSums(xyCors) * colSums(Y_scaled^2 * t(xRowSum2))
        
        # Calc for star 4
        star4 <- rowSums(xRowSum2^2)
        
        # Calc for dagger 1
        dagger1 <- rowSums(xyCors) * crossprod(Y_scaled^3, xRowSum)
        
        # Calc for dagger 2
        dagger2 <- colSums(xRowSum * t(xRowSum2) * Y_scaled)
      }
      
    }
    
    allvars <- (star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / 
      (n - 1)
    corsums <- as.vector(rowSums(xyCors))
    zstats <- sqrt(n) * corsums / sqrt(allvars)
    pvals <- pnorm(zstats, lower.tail = FALSE)
    
    if (justpvals)
      return(pvals)
    
    successes <- bh_reject(pvals, alpha)
    
    # Return indices
    
    if (test_X) {
      Anew <- Xindx[successes]
    } else {
      Anew <- Yindx[successes]
    }
    
    return(Anew)
    
  }
  
  initialize3 <- function(u){
    if (u <= dx) {
      cor_to_u <- cross_cor(X_indx=u)
      fischer_tranformed_cor <- atanh(cor_to_u)*sqrt(n-3)
      pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
      successes <- Yindx[bh_reject(pvals, alpha)]
    } else {
      if (u > dx) {
        cor_to_u <- cross_cor(Y_indx=u-dx)
        fischer_tranformed_cor <- atanh(cor_to_u)*sqrt(n-3)
        pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
        successes <- Xindx[bh_reject(pvals, alpha)]
      } 
    }
    return(successes)
  }
  
  initialize <- function (...) {
    
    if (initializeMethod == 3) {
      return(initialize3(...))
    } else {
      stop('only initializeMethod = 3 supported\n')
    }
    
    
  }
  
  update <- function (...) {
    
    if (updateMethod == 5)
      return(update5(...))
    
    if (!updateMethod %in% c(5, 7)) {
      stop('only updateMethod = 5 or 7 supported\n')
    }
    
  }
  
  #-------------------------------------------------------------------------------
  # Extractions set-up
  
  cat("Beginning method.\n\n")
  
  # Getting node orders
  Ysum <- Y_scaled %*% rep(1,dy)
  Xsum <- X_scaled %*% rep(1,dx)
  cor_X_to_Ysums <- as.vector(t(Ysum) %*% X_scaled)
  cor_Y_to_Xsums <- as.vector(t(Xsum) %*% Y_scaled)
  
  remainingY <- Yindx[order(cor_Y_to_Xsums, decreasing = TRUE)]
  remainingX <- Xindx[order(cor_X_to_Ysums, decreasing = TRUE)]
  
  # Initializing control variables
  didX <- TRUE
  loopCount <- 0
  OL_count <- 0
  Dud_count <- 0
  comms <- rep(list(integer(0)), dx + dy)
  initial.sets <- comms
  final.sets <- comms
  did_it_cycle <- logical(dx + dy)
  update_info <- comms
  nits <- numeric(dx + dy)
  plugged <- logical(dx + dy)
  
  # Making save string
  starttime <- paste0(unlist(strsplit(as.character(starttime), " ", fixed = TRUE)), collapse = "_")
  starttime <- gsub(":", "-", starttime)
  
  if (!is.null(saveDir)) {
    if (!is.null(tag)) {
      fn <- file.path(saveDir, paste0(tag, ".RData"))
    } else {
      fn <- file.path(saveDir, paste0("unnamed", starttime, ".RData"))
    }
  }
  
  if (!is.null(saveDir)) {
    cat("doing test save\n")
    cat("fn is", fn, "\n")
    # Test save
    testsave <- "testsave"
    save(testsave, file = fn)
  }
  
  rm(comms)
  
  while (length(c(remainingX, remainingY)) > 0) {
    
    if (var(c(length(final.sets), length(update_info), length(initial.sets))) != 0)
      break
    
    cat("\n#-----------------------------------------\n\n")
    loopCount <- loopCount + 1
    cat("loopCount", loopCount, "\n\n")
    cat("length remainingX = ", length(remainingX), "\n")
    cat("length remainingY = ", length(remainingY), "\n")
    cat("OL_count = ", OL_count, "\n")
    cat("Dud_count = ", Dud_count, "\n")
    
    # Initializing
    if (length(remainingX) == 0 || didX & length(remainingY) > 0) {
      comm_indx <- remainingY[1]
      remainingY <- remainingY[-1]
      didX <- FALSE
      comm_indxY <- which(Yindx == comm_indx)
      cat("comm_starter = ", comm_indx, "(Y =", paste0(comm_indxY, ")"), "\n")
    } else {
      comm_indx <- remainingX[1]
      remainingX <- remainingX[-1]
      didX <- TRUE
      comm_indxX <- which(Xindx == comm_indx)
      cat("comm_starter = ", comm_indx, "(X =", paste0(comm_indxX, ")"), "\n")
    }
    
    
    # Update 0
    cat("getting, checking initial sets\n")
    B0x <- initialize(comm_indx)
    if (length(B0x) <= 1) {
      Dud_count <- Dud_count + 1
      
      if (!is.null(saveDir)) {
        # Saving status
        save(loopCount, initial.sets, final.sets, update_info,
             OL_count, Dud_count,
             file = fn)
      }
      
      if (Dud_count > Dud_tol)
        break
      next
    }
    if (updateMethod == 7) {
      B0y <- update5(B0x, comm_indx)
    } else {
      B0y <- update(B0x, comm_indx)
    }
    
    # Removing B0 from remaining
    if (throwInitial) {
      remainingX <- setdiff(remainingX, c(B0x, B0y))
      remainingY <- setdiff(remainingY, c(B0x, B0y))
    }
    
    # Initializing update loop
    B_oldx <- B0x
    B_oldy <- B0y
    B_old <- c(B_oldx, B_oldy)
    initial.sets[[comm_indx]] <- B_old
    B_new <- c(Xindx, Yindx)
    chain <- list(B_old)
    consec_jaccards <- NULL
    consec_sizes <- list(c(length(B_oldx), length(B_oldy)))
    if(is.function(bmd_index)){
      consec_composition <- list(table(sapply(B_old,bmd_index)))
    } else {
      consec_composition <- NULL
    }
    
    found_cycle <- found_break <- NULL
    mean_jaccards <- NULL ## add all these to update_info
    itCount <- 0
    cycledSets <- NULL
    if (updateOutput) {
      cat(paste0("Initial set ",
                 "is size ", length(B_old), 
                 " (", length(B_oldx), ", ", length(B_oldy), ")\n\n", sep=""))
    }
    
    repeat {
      
      itCount <- itCount + 1
      
      if (updateMethod != 7) {
        B_newx <- update(B_oldy, B_oldx)
        if (length(B_newx) >= 1) {
          B_newy <- update(B_newx, B_oldy)
        } else {
          B_newy <- integer(0)
          break
        }
        B_new <- c(B_newx, B_newy)
        
        if (length(B_newy) == 0)
          break
      } else {
        
        if (comm_indx > dx) {
          Xpvals <- update5(B_oldy, justpvals = TRUE)
          Ypvals <- update5(B_oldx, justpvals = TRUE)
          B_new <- bh_reject(c(Xpvals, Ypvals), alpha)
          B_newx <- B_new[B_new <= dx]
          B_newy <- B_new[B_new > dx]
        } else {
          Xpvals <- update5(B_oldx, justpvals = TRUE)
          Ypvals <- update5(B_oldy, justpvals = TRUE)
          B_new <- bh_reject(c(Xpvals, Ypvals), alpha)
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
      if(is.function(bmd_index)){
        consec_composition <- c(consec_composition,
                                list(table(sapply(B_new, bmd_index))))
      }
      if (updateOutput) {
        cat(paste0("Update ", itCount, 
                   " is size ", length(B_new), 
                   " (", length(B_newx), ", ", length(B_newy), "), ",
                   "jaccard to last is ", round(consec_jaccard, 3), ", ",
                   "mean jaccard along chain is ", round(mean(jaccards), 3), "\n", sep=""))
      }
      
      # Checking for cycles (4.4.1 in CCME paper)
      if (jaccard(B_new, B_old) > 0) { # Otherwise loop will end naturally
        
        if (sum(jaccards == 0) > 0) { # Cycle has been found
          
          found_cycle[itCount] <- TRUE
          did_it_cycle[comm_indx] <- TRUE
          
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
      
    } # From Updates
    
    nits[comm_indx] <- itCount
    
    if (length(B_newx) * length(B_newy) == 0)
      B_new = integer(0)
    
    #comms[[comm_indx]] <- B_new
    final.sets[[comm_indx]] <- B_new
    plugged[comm_indx] <- length(B_new) > 0
    if (sum(unlist(lapply(final.sets[plugged], length)) == 0) > 0)
      break
    update_info[[comm_indx]] <- list("mean_jaccards" = mean_jaccards,
                                     "consec_jaccards" = consec_jaccards,
                                     "consec_sizes"= consec_sizes,
                                     "consec_composition"=consec_composition,
                                     "found_cycle" = found_cycle,
                                     "found_break" = found_break)
    remainingX <- setdiff(remainingX, B_new)
    remainingY <- setdiff(remainingY, B_new)
    
    if (sum(plugged) > 0) {
      
      # Checking overlap with previous sets
      cat("\nchecking overlap...\n\n")
      OL_check <- unlist(lapply(final.sets[plugged], function (C) jaccard(B_new, C)))
      #OL_check[is.na(OL_check)] <- 1
      OL_check[which(plugged) == comm_indx] <- 1
      if (sum(OL_check < 1 - OL_thres) > 0)
        OL_count <- OL_count + 1
      
    }
    
    comm_nodes <- unique(unlist(final.sets))
    cat(paste0(sum(comm_nodes <= dx), " X vertices in communities.\n"))
    cat(paste0(sum(comm_nodes > dx), " Y vertices in communities.\n"))
    
    if (!is.null(saveDir)) {
      # Saving status
      save(loopCount, initial.sets, final.sets, update_info,
           OL_count, Dud_count,
           file = fn)
    }
    
    current_time <- proc.time()[3] - start_second
    if (current_time > time_limit)
      break
    
    if (OL_count > OL_tol)
      break
    
  }
  
  #-----------------------------------------------------------------------------
  #  Clean-up and return -------------------------------------------------------
  
  endtime <- Sys.time()
  report <- list(loopCount, OL_count, Dud_count, starttime, endtime)
  names(report) <- c("loopCount", "OL_count", "Dud_count", "starttime", "endtime")
  
  # Removing blanks and trivial sets
  cat("removing overlap and unpacking comms...\n")
  nonNullIndxs <- which(unlist(lapply(final.sets, length)) > 0)
  if (length(nonNullIndxs) == 0) {
    
    returnList <- list("communities" = NULL,
                       "background" = list("X_bg" = 1:dx,
                                           "Y_bg" = (dx + 1):(dx + dy)),
                       "commzs" = NULL,
                       "initial.sets" = initial.sets, 
                       "final.sets" = final.sets,
                       "nits" = nits,
                       "report" = report,
                       "communities_before_OLfilt" = NULL,
                       "OLfilt" = NULL,
                       "did_it_cycle" = NULL,
                       "update_info" = update_info,
                       "finalIndxs" = NULL,
                       "nonNullIndxs" = nonNullIndxs)
    return(returnList)
    
    
  }
  nonNullComms <- final.sets[nonNullIndxs]
  OLfilt <- filter_overlap(nonNullComms, tau = OL_thres, inv.length = inv.length)
  finalComms <- OLfilt$final_comms
  finalIndxs <- nonNullIndxs[OLfilt$kept_comms]
  final_did_it_cycle <- did_it_cycle[finalIndxs]
  
  X_sets <- lapply(finalComms, function (C) C[C <= dx])
  Y_sets <- lapply(finalComms, function (C) C[C > dx])
  X_bg <- setdiff(Xindx, unlist(X_sets))
  Y_bg <- setdiff(Yindx, unlist(Y_sets))
  final_ncomms <- length(X_sets)
  
  returnList <- list("communities" = list("X_sets" = X_sets,
                                          "Y_sets" = Y_sets),
                     "background" = list("X_bg" = X_bg,
                                         "Y_bg" = Y_bg),
                     "initial.sets" = initial.sets, 
                     "final.sets" = final.sets,
                     "nits" = nits,
                     "report" = report,
                     "communities_before_OLfilt" = nonNullComms,
                     "OLfilt" = OLfilt,
                     "did_it_cycle" = final_did_it_cycle,
                     "update_info" = update_info,
                     "finalIndxs" = finalIndxs,
                     "nonNullIndxs" = nonNullIndxs)
  return(returnList)
  
}
