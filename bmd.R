source("makeVars.R")
source("stdize.R")
bmd <- function (X, Y, alpha = 0.05, OL_thres = 0.9, tag = NULL, saveDir = getwd(),
                 updateOutput = TRUE, throwInitial = TRUE, OL_tol = Inf, Dud_tol = Inf, time_limit = 18000,
                 updateMethod = 2, initializeMethod = 2, inv.length = TRUE, add_rate = 1, return_zs = TRUE) {

  if (FALSE) {
    alpha = 0.05
    OL_thres = 0.9
    tag = NULL
    saveDir = getwd()
    OL_tol = Inf
    Dud_tol = Inf
    time_limit = 18000
    updateMethod = 2
    initializeMethod = 2
    updateOutput = TRUE
    throwInitial = TRUE
    inv.length = TRUE
    time_limit = Inf
    add_rate = 1
    return_zs = TRUE
  }
  
  start_second <- proc.time()[3]
  starttime <- Sys.time()
  if (!dir.exists(saveDir))
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
  
  bhy <-
    function(pvals, alpha = 0.05){
      
      # Sorted p-vals
      sp = sort(pvals)
      
      # Save original order of p-vals
      ord = order(pvals)
      
      # Find bhy cutoff
      nums = 1:length(pvals)
      cms = cumsum(1/nums)
      
      # Find which p-vals are less than bh cutoff
      under = sp < (nums/(length(pvals)*cms)*alpha)
      
      # Return indices of significant p-vals
      if(sum(under) == 0){
        return(c())
      }else{
        cutoff = max(which(under))
        return(ord[1:cutoff])
      }
    }
  
  bh_reject <- function (pvals, alpha, conserv = TRUE) {
    
    m <- length(pvals)
    
    if (!conserv) {
      pvals_adj <- m * pvals / rank(pvals)
    } else {
      mults <- cumsum(1 / c(1:m))
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
  
  Xindx <- 1:dx
  Yindx <- (dx + 1):(dx + dy)
  
  # Defining pval function
  pvalFun <- function (B) {
    
    side2 <- min(B) > dx
    
    if (side2) {
      
      B_match <- match(B, Yindx)
      R_mat <- crossprod(Y_scaled[ , B_match], X_scaled) / (n - 1)
      stats <- colSums(R_mat)
      if (length(B) > 1) {
        EmpCov <- cor(Y[ , B_match])
      } else {
        EmpCov <- 1
      }
      varSum <- sum(EmpCov)
      
    } else { 
      
      B_match <- match(B, Xindx)
      R_mat <- crossprod(X_scaled[ , B_match], Y_scaled) / (n - 1)
      stats <- colSums(R_mat)
      if (length(B) > 1) {
        EmpCov <- cor(X[ , B_match])
      } else {
        EmpCov <- 1
      }
      varSum <- sum(EmpCov)
      
    }
    
    pvals <- pnorm(sqrt(n) * stats, sd = sqrt(varSum), lower.tail = FALSE)
    
    return(pvals)
    
  }
  
  update1 <- function (B, B0 = NULL) {
    
    side2 <- min(B) > dx
      
    pvals <- pvalFun(B)
    
    if (side2) {
      B_new <- bh_reject(pvals, alpha)
    } else {
      B_new <- Yindx[bh_reject(pvals, alpha)]
    }
    
    return(B_new)
    
  }
  
  scoreCalcs <- function (B, A) {
    
    I <- length(A)
    J <- length(B)
    
    if (max(A) <=  dx) {
      RA <- sum(cor(as.matrix(X[ , A])))
      RB <- sum(cor(as.matrix(Y[ , match(B, Yindx)])))
      Stat <- n^(1 / 2) * sum(cor(as.matrix(X[ , A]), as.matrix(Y[ , match(B, Yindx)])))
    } else { 
      RA <- sum(cor(as.matrix(Y[ , match(A, Yindx)])))
      RB <- sum(cor(as.matrix(X[ , B])))
      Stat <- n^(1 / 2) * sum(cor(as.matrix(Y[ , match(A, Yindx)]), as.matrix(X[ , B])))
    }
    
    VAB <- J * RA + I * RB - I * J
      
    return(c(RA, RB, VAB, Stat))
    
  }
  
  update2 <- function (B, A = NULL) {
    
    test_X <- min(B) > dx
    nFixd <- length(B)
    
    chunkSize <- 500
    
    
    if (test_X) {
      
      # Getting fixd matrix and calcs
      fixdIndx <- match(B, Yindx)
      fixdMat <- as.matrix(Y[ , fixdIndx])
      
      # Getting test stats and indices
      all_stats <- n^(1 / 2) * colSums(cor(fixdMat, X))
      nTest <- dx
      testOrder <- order(all_stats, decreasing = TRUE)
      
      
    } else {
      
      # Getting fixd matrix and calcs
      fixdIndx <- match(B, Xindx)
      fixdMat <- as.matrix(X[ , fixdIndx])
      
      # Getting test stats and indices
      all_stats <- n^(1 / 2) * colSums(cor(fixdMat, Y))
      nTest <- dy
      testOrder <- order(all_stats, decreasing = TRUE)
      
    }
    
    nChunks <- ceiling(nTest / chunkSize)
    c_pos <- 1
    total_includes <- numeric(nTest)
    total_zs <- numeric(nTest)
    last_z <- 0
    
    for (c in 1:nChunks) {
      
      cat("........chunk", c, "\n")
      
      indxs <- c_pos:min(nTest, (c_pos + chunkSize - 1))
      if (test_X) {
        A1 <- Xindx[testOrder[indxs]]
      } else {
        A1 <- Yindx[testOrder[indxs]]
      }
      zs <- numeric(length(indxs))
      vars <- zs
      stats <- zs
      
      for (t in seq_along(A1)) {
        
        score_comps <- scoreCalcs(B, A1[1:t])
        stats[t] <- score_comps[4]
        vars[t] <- score_comps[3]
        zs[t] <- stats[t] / sqrt(vars[t])
        
      }
      
      includes <- sign(diff(c(last_z, zs)))
      total_includes[indxs] <- includes
      
      if (sum(includes) != length(includes))
        break
      
      c_pos <- c_pos + chunkSize
      last_z <- tail(zs, 1)
      
    }
    
    last_in <- min(which(total_includes == -1)) - 1
    
    if (test_X) {
      Anew <- Xindx[testOrder[1:last_in]]
    } else {
      Anew <- Yindx[testOrder[1:last_in]]
    }
    
    return(Anew)
    
  }
  
  update2fast <- function(B, A = NULL, final = FALSE, num_new = 0) {
    
    test_X <- min(B) > dx
    J <- length(B)
    
    chunkSize <- 60
    
    
    if (test_X) {
      
      # Getting fixd matrix and calcs
      fixdIndx <- match(B, Yindx)
      fixdMat <- as.matrix(Y[ , fixdIndx])
      
      # Getting test stats and indices
      all_stats <- n^(1 / 2) * colSums(cor(fixdMat, X))
      nTest <- dx
      testOrder <- order(all_stats, decreasing = TRUE)
      
      
    } else {
      
      # Getting fixd matrix and calcs
      fixdIndx <- match(B, Xindx)
      fixdMat <- as.matrix(X[ , fixdIndx])
      
      # Getting test stats and indices
      all_stats <- n^(1 / 2) * colSums(cor(fixdMat, Y))
      nTest <- dy
      testOrder <- order(all_stats, decreasing = TRUE)
      
    }
    
    # Fixed quantities
    nChunks <- ceiling(nTest / chunkSize)
    c_break <- nChunks
    total_includes <- logical(nTest)
    total_zs <- numeric(nTest)
    RB <- sum(cor(fixdMat))
    
    # Rolling quantities
    c_pos <- 1
    lastV <- 0
    Rhat_A0_B <- 0
    RA0 <- 0
    A0 <- integer(0)
    do_more <- FALSE
    
    for (c in 1:nChunks) {
      
      cat("........chunk", c, "\n")

      indxs <- c_pos:min(nTest, c_pos + chunkSize - 1)
      
      A1 <- testOrder[indxs]
      K1 <- length(A1)
      K0 <- length(A0)
      
      R_k_A0 <- 0
      
      if (test_X) {
        if (K0 > 0) {
          R_k_A0 <- rowSums(cor(X[ , A1], X[ , A0]))
        }
        Rhat_k_B <- rowSums(cor(X[ , A1], Y[ , fixdIndx]))
        RA1 <- cor(X[ , A1])
      } else {
        if (K0 > 0) {
          R_k_A0 <- rowSums(cor(Y[ , A1], Y[ , A0]))
        }
        Rhat_k_B <- rowSums(cor(Y[ , A1], X[ , fixdIndx]))
        RA1 <- cor(Y[ , A1])
      }
      
      RA1[lower.tri(RA1)] <- 0; diag(RA1) <- 0
      R_k_A1_back <- colSums(RA1)
      A01_vars <- 2 * (R_k_A0 + R_k_A1_back) + 1
      RA0k <- RA0 + cumsum(A01_vars)
      V_A0k_B <- J * RA0k + indxs * RB - indxs * J
      Rhat_A0k_B <- Rhat_A0_B + cumsum(Rhat_k_B)
      
      Vratios <- V_A0k_B / c(lastV, V_A0k_B[-K1])
      Rsums <- c(Rhat_A0_B, Rhat_A0k_B[-K1])
      
      zs <- (Rhat_k_B + Rsums) / sqrt(V_A0k_B)
      total_zs[indxs] <- zs
      if (c == c_break)
        break
      
      if (!do_more) {
        includes <- Rhat_k_B > Rsums * (sqrt(Vratios) - 1)
        if (c == 1) {
          includes[1] <- TRUE
        }
        total_includes[indxs] <- includes
      }
      
      if (!do_more && sum(includes) != length(includes)) {
        if (!final) {
          break
        } else {
          do_more <- TRUE
          final_spot <- min(which(includes == FALSE)) - 1
          left_in_chunk <- chunkSize - final_spot
          needed <- num_new - left_in_chunk
          if (needed <= 0) {
            break
          } else {
            if (needed + c_pos + chunkSize - 1 > nTest)
              needed <- nTest - c_pos - chunkSize + 1
            chunks_needed <- ceiling(needed / chunkSize)
            c_break <- c + chunks_needed
          }
        }
      }
        
      
      c_pos <- c_pos + chunkSize
      lastV <- V_A0k_B[K1]
      Rhat_A0_B <- Rhat_A0k_B[K1]
      RA0 <- RA0k[K1]
      A0 <- c(A0, A1)
      
    }
    
    last_in <- min(which(total_includes == FALSE)) - 1
    if (final)
      last_in2 <- min(min(which(total_zs == 0)) - 1, 1)
      
    if (test_X) {
      Anew <- Xindx[testOrder[1:last_in]]
      if (final)
        Anew2 <- Xindx[testOrder[1:last_in2]]
    } else {
      Anew <- Yindx[testOrder[1:last_in]]
      if (final)
        Anew2 <- Yindx[testOrder[1:last_in2]]
    }
    
    if (final) {
      return(list("updated" = Anew,
                  "updated2" = Anew2,
                  "zs" = total_zs[1:length(Anew2)]))
    } else {
      return(Anew)
    }
    
  }
  
  update4 <- function (B, A) {
    
    test_X <- max(A) <= dx
    nTest <- length(A)
    nFixd <- length(B)
    
    if (nFixd == 0)
      return(integer(0))
    
    if (test_X) {
      
      # Getting indices
      testIndx <- match(A, Xindx)
      fixdIndx <- match(B, Yindx)
      set_ind <- Xindx %in% A
      
      # Getting matrices
      testMat <- as.matrix(X[ , testIndx])
      fixdMat <- as.matrix(Y[ , fixdIndx])
      
    } else {
      
      # Getting indices
      testIndx <- match(A, Yindx)
      fixdIndx <- match(B, Xindx)
      set_ind <- Yindx %in% A
      
      # Getting matrices
      testMat <- as.matrix(Y[ , testIndx])
      fixdMat <- as.matrix(X[ , fixdIndx])
      
    }
    
    # Getting variance sums
    RB <- sum(cor(fixdMat))
    RA <- sum(cor(testMat))
    
    # Getting big sum
    Rhat_A_B <- n^(1 / 2) * sum(cor(fixdMat, testMat))
    Rhat_A_B_i <- Rhat_A_B - set_ind
    
    # Getting all sum stats and var means
    if (test_X) {
      Rhat_i_B <- n^(1 / 2) * colSums(cor(fixdMat, X))
      R_i_A <- colSums(cor(testMat, X)) - set_ind
    } else {
      Rhat_i_B <- n^(1 / 2) * colSums(cor(fixdMat, Y))
      R_i_A <- colSums(cor(testMat, Y)) - set_ind
    }
    RA_i <- RA - set_ind * (2 * R_i_A + 1)
    
    Fi <- 1 + (nFixd * (2 * R_i_A) + RB) / (nFixd * RA_i + (nTest - set_ind) * (RB - nFixd))
    successes <- Rhat_i_B >= Rhat_A_B_i * (sqrt(Fi) - 1)
    score_ratios <- Rhat_i_B / (Rhat_A_B_i * (sqrt(Fi) - 1))
    
    
    # Return indices
    
    if (test_X) {
      Anew <- Xindx[successes]
      testIndx <- Xindx[testIndx]
    } else {
      Anew <- Yindx[successes]
      testIndx <- Yindx[testIndx]
    }
    
    if (length(testIndx) == 1 && length(Anew) > 0)
      Anew <- c(testIndx, Anew)
    
    return(Anew)
    
  }
  
  update5 <- function (B, A = NULL, bhy = FALSE) {
    
    if (length(B) == 0)
      return(integer(0))
    
    test_X <- min(B) > dx
    nFixd <- length(B)
    
    if (test_X) {
      
      # Getting fixed matrix
      fixdIndx <- match(B, Yindx)
      fixdMat <- as.matrix(Y_scaled[ , fixdIndx])

      # Calculating the variances
      {
        # General calcs
        xyCors <- cor(X_scaled, fixdMat)
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
      fixdMat <- as.matrix(X_scaled[ , fixdIndx])
      
      # Calculating the variances
      {
        # General calcs
        xyCors <- cor(Y_scaled, fixdMat)
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
    
    if (bhy) {
      successes <- bhy(pvals, alpha)
    } else {
      successes <- bh_reject(pvals, alpha)
    }
    
    
    # Return indices
    
    if (test_X) {
      Anew <- Xindx[successes]
    } else {
      Anew <- Yindx[successes]
    }
    
    return(Anew)
    
  }
  
  update6 <- function (B, A = NULL) {
    
    if (length(B) == 0)
      return(integer(0))
    
    test_X <- min(B) > dx
    nFixd <- length(B)
    
    if (test_X) {
      
      # Getting fixed matrix
      fixdIndx <- match(B, Yindx)
      fixdMat <- as.matrix(Y[ , fixdIndx])
      Uis <- X
      cormeans <- as.vector(rowMeans(cor(X, fixdMat)))
      
    } else {
      
      # Getting fixed matrix
      fixdIndx <- match(B, Xindx)
      fixdMat <- as.matrix(X[ , fixdIndx])
      Uis <- Y
      cormeans <- as.vector(rowMeans(cor(Y, fixdMat)))
      
    }
    
    Uis <- stdize(t(Uis))
    M_A <- stdize(t(fixdMat))
    
    
    
    allvars <- makeVars(Uis, M_A)
    zstats <- cormeans / sqrt(allvars)
    pvals <- pnorm(zstats, lower.tail = FALSE)
    successes <- bhy(pvals, alpha)
    
    
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
      cor_to_u <- cor(X_scaled[,u], Y_scaled)
      fischer_tranformed_cor <- atanh(cor_to_u)*sqrt(n-3)
      pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
      successes <- Yindx[bh_reject(pvals, alpha)]
    } else {
      if (u > dx) {
        cor_to_u <- cor(Y_scaled[,u-dx], X_scaled)
        fischer_tranformed_cor <- atanh(cor_to_u)*sqrt(n-3)
        pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
        successes <- Xindx[bh_reject(pvals, alpha)]
      } 
    }
    return(successes)
  }
  
  initialize <- function (...) {
    
    if (initializeMethod == 1)
      return(update1(...))
    
    if (initializeMethod == 2)
      return(update2fast(...))
    
    if (initializeMethod == 3)
      return(initialize3(...))
    

  }
  
  update <- function (...) {
    
    if (updateMethod == 1)
      return(update1(...))
    
    if (updateMethod == 2)
      return(update2fast(...))
    
    if (updateMethod == 4)
      return(update4(...))
    
    if (updateMethod == 5)
      return(update5(...))
    
    if (updateMethod == 6)
      return(update6(...))
      
  }
  
  #-------------------------------------------------------------------------------
  
  # Extractions set-up

  cat("Beginning method.\n\n")
  
  # Getting node orders
  Yvar <- apply(Y, 2, var)
  remainingY <- Yindx[order(Yvar, decreasing = TRUE)]
  Xvar <- apply(X, 2, var)
  remainingX <- Xindx[order(Xvar, decreasing = TRUE)]

  # Initializing control variables
  didX <- TRUE
  loopCount <- 0
  OL_count <- 0
  Dud_count <- 0
  comms <- rep(list(integer(0)), dx + dy)
  initial.sets <- comms
  final.sets <- comms
  did_it_cycle <- logical(length(comms))
  update_info <- comms
  nits <- numeric(dy)
  
  # Making save string
  starttime <- paste0(unlist(strsplit(as.character(starttime), " ", fixed = TRUE)), collapse = "_")
  starttime <- gsub(":", "-", starttime)
  if (!is.null(tag)) {
    fn <- file.path(saveDir, paste0(tag, ".RData"))
  } else {
    fn <- file.path(saveDir, paste0("unnamed", starttime, ".RData"))
  }
  
  cat("doing test save\n")
  cat("fn is", fn, "\n")
  # Test save
  testsave <- "testsave"
  save(testsave, file = fn)
  
  while (length(c(remainingX, remainingY)) > 0) {
    
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
      
      # Saving status
      save(loopCount, comms, initial.sets, final.sets, 
           OL_count, Dud_count,
           file = fn)
      
      if (Dud_count > Dud_tol)
        break
      next
    }
    B0y <- update(B0x, comm_indx)
    
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
    consec_sizes <- list()
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
      
      #if (itCount == 50)
      #  break
      
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
      
      consec_jaccard <- jaccard(B_new, B_old)
      jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
      found_cycle <- c(found_cycle, FALSE)
      found_break <- c(found_break, FALSE)
      consec_jaccards <- c(consec_jaccards, consec_jaccard)
      mean_jaccards <- c(mean_jaccards, mean(jaccards))
      consec_sizes <- c(consec_sizes, list(c(length(B_newx), length(B_newy))))
      
      if (updateOutput) {
        cat(paste0("Update ", itCount, 
                   " is size ", length(B_new), 
                   " (", length(B_newx), ", ", length(B_newy), "), ",
                   "jaccard to last is ", round(consec_jaccard, 3), ", ",
                   "mean jaccard along chain is ", round(mean(jaccards), 3), "\n", sep=""))
      }
      
      #TODO: store all the "last" jaccards in a vector to be returned
      
      
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
            B_new <- NULL
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
    
    comms[[comm_indx]] <- B_new
    final.sets[[comm_indx]] <- B_new
    update_info[[comm_indx]] <- list("mean_jaccards" = mean_jaccards,
                                     "consec_jaccards" = consec_jaccards,
                                     "consec_sizes"= consec_sizes,
                                     "found_cycle" = found_cycle,
                                     "found_break" = found_break)
    remainingX <- setdiff(remainingX, B_new)
    remainingY <- setdiff(remainingY, B_new)
    
    # Checking overlap with previous sets
    cat("\nchecking overlap...\n\n")
    OL_check <- unlist(lapply(comms, function (C) jaccard(B_new, C)))
    OL_check[is.na(OL_check)] <- 1
    OL_check[comm_indx] <- 1
    if (sum(OL_check < 1 - OL_thres) > 0)
      OL_count <- OL_count + 1
    
    comm_nodes <- unique(unlist(comms))
    cat(paste0(sum(comm_nodes <= dx), " X vertices in communities.\n"))
    cat(paste0(sum(comm_nodes > dx), " Y vertices in communities.\n"))
    
    # Saving status
    save(loopCount, chain, comms, initial.sets, final.sets, 
         OL_count, Dud_count,
         file = fn)
    
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
  nonNullIndxs <- which(unlist(lapply(comms, length)) > 0)
  if (length(nonNullIndxs) == 0) {
    
    returnList <- list("communities" = comms,
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
  nonNullComms <- comms[nonNullIndxs]
  OLfilt <- filter_overlap(nonNullComms, tau = OL_thres, inv.length = inv.length)
  finalComms <- OLfilt$final_comms
  finalIndxs <- nonNullIndxs[OLfilt$kept_comms]
  final_did_it_cycle <- did_it_cycle[finalIndxs]

  X_sets <- lapply(finalComms, function (C) C[C <= dx])
  Y_sets <- lapply(finalComms, function (C) C[C > dx])
  X_bg <- setdiff(Xindx, unlist(X_sets))
  Y_bg <- setdiff(Yindx, unlist(Y_sets))
  final_ncomms <- length(X_sets)
  X_zs <- rep(list(NULL), final_ncomms)
  Y_zs <- rep(list(NULL), final_ncomms)
  
  if (return_zs) {
    
    for (c in 1:final_ncomms) {
      
      cat("Calculating final zs for comm", c, "\n")
      
      XL <- length(X_sets[[c]])
      YL <- length(Y_sets[[c]])
      num_newX <- max(5, ceiling(add_rate * XL))
      num_newY <- max(5, ceiling(add_rate * YL))
      updateObjX <- update2fast(Y_sets[[c]], final = TRUE, num_new = num_newX)
      updateObjY <- update2fast(X_sets[[c]], final = TRUE, num_new = num_newX)
      X_zs[[c]] <- updateObjX[2:3]
      Y_zs[[c]] <- updateObjY[2:3]
      
    }
  }
      
  
  
  
  returnList <- list("communities" = list("X_sets" = X_sets,
                                          "Y_sets" = Y_sets),
                     "background" = list("X_bg" = X_bg,
                                         "Y_bg" = Y_bg),
                     "commzs" = list("X_zs" = X_zs,
                                     "Y_zs" = Y_zs),
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