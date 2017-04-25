extract <- function (indx, interact = FALSE, print_output = verbose) {
  
  if (print_output) {
    cat("\n#-----------------------------------------\n\n")
    cat("trying indx", indx, "\n")
  }
  
  # If you want to interact with trackers, get the current tracking info
  if (interact) {
    
    # Seeing if stop has been triggered
    stop_extracting <- as.logical(readLines(stop_fn))
    if (stop_extracting) return(list(report = "stop_extracting"))
    
    # Seeing if indx already in a comm
    comm_files <- list.files(comm_dn, full.names = TRUE)
    if (length(comm_files) > 0) {
      comms <- lapply(comm_files, readLines)
    } else {
      comms <- numeric(0)
    }
    clustered <- as.numeric(unique(unlist(comms)))
    if (indx %in% clustered) return(list(report = "indx_clustered"))
    
    OL_count <- as.numeric(readLines(OL_fn))
    Dud_count <- as.numeric(readLines(Dud_fn))
    
  }
  
  # Print some output
  if (print_output && interact) {
    cat("extraction", which(extractord == indx), "of", length(extractord), "\n\n")
    cat("OL_count = ", OL_count, "\n")
    cat("Dud_count = ", Dud_count, "\n")
    cat(paste0(sum(clustered <= dx), " X vertices in communities.\n"))
    cat(paste0(sum(clustered > dx), " Y vertices in communities.\n"))
  }
  
  # Get general B01 (regardless of indx)
  if (Cpp) {
    B01 <- initialize1Cpp(indx)
  } else {
    B01 <- initialize1R(indx)
  }
  
  if (length(B01) > 1) {
    
    # Convert indices if necessary
    if (indx <= dx) B01 <- Yindx[B01]
    
    # Half-update
    if (Cpp) {
      pvals2 <- pvalsCpp(B01)
      B02 <- bh_rejectC(pvals2, alpha, conserv = TRUE)
    } else {
      pvals2 <- pvalsR(B01)
      B02 <- bh_rejectR(pvals2, alpha, conserv = TRUE)
    }
    
  } else {
    B02 <- integer(0)
  }
  
  # Convert indices if necessary
  if (indx > dx) B02 <- Yindx[B02]
  
  # Check for dud
  if (length(B01) * length(B02) <= 1) {
    if (interact) {
      Dud_count <- Dud_count + 1
      writeLines(as.character(Dud_count), Dud_fn)
    }
    return(NULL)
  }
  
  # Assign B0x and B0y according to indx
  if (indx > dx) {
    B0x <- B01; B0y <- B02
  } else {
    B0x <- B02; B0y <- B01
  }
  
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
  
  # Extraction loop
  repeat {
    
    itCount <- itCount + 1
    
    if (updateMethod == 1) {
      
      pvalsXY <- if (Cpp) {
        c(pvalsCpp(B_oldy), pvalsCpp(B_oldx))
      } else {
        c(pvalsR(B_oldy), pvalsR(B_oldx))
      }
      B_new <- if (Cpp) {
        bh_rejectC(pvalsXY, alpha, conserv = TRUE)
      } else {
        bh_rejectR(pvalsXY, alpha)
      }
      B_newx <- B_new[B_new <= dx]; B_newy <- B_new[B_new > dx]
      
    }
    
    if (updateMethod == 2) {
      
      if (indx > dx) {
        
        # Update X nodes first
        pvalsx <- if (Cpp) pvalsCpp(B_oldy) else pvalsR(B_oldy)
        B_newx <- if(Cpp) {
          bh_rejectC(pvalsx, alpha, conserv = TRUE)
        } else {
          bh_rejectR(pvalsx, alpha, conserv = TRUE)
        }
        if (length(B_newx) <= 1) {
          B_newy <- integer(0)
          break
        }
        pvalsy <- if (Cpp) pvalsCpp(B_newx) else pvalsR(B_newx)
        B_newy <- if (Cpp) {
          Yindx[bh_rejectC(pvalsy, alpha, conserv = TRUE)]
        } else {
          Yindx[bh_rejectR(pvalsy, alpha, conserv = TRUE)]
        }
        
      } else {
        
        # Update Y nodes first
        pvalsy <- if (Cpp) pvalsCpp(B_oldx) else pvalsR(B_oldx)
        B_newy <- if (Cpp) {
          Yindx[bh_rejectC(pvalsy, alpha, conserv = TRUE)]
        } else {
          Yindx[bh_rejectR(pvalsy, alpha, conserv = TRUE)]
        }
        if (length(B_newy) <= 1) {
          B_newx <- integer(0)
          break
        }
        pvalsx <- if (Cpp) pvalsCpp(B_newy) else pvalsR(B_newy)
        B_newx <- if (Cpp) {
          bh_rejectC(pvalsx, alpha, conserv = TRUE)
        } else {
          bh_rejectR(pvalsx, alpha, conserv = TRUE)
        }
        
      }
        
      B_new <- c(B_newx, B_newy)
      
    }
    
      
    if (length(B_newy) * length(B_newx) == 0) {
      B_newy <- B_newx <- integer(0)
      break
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
  if (length(B_newx) * length(B_newy) == 0) {
    B_new <- integer(0)
    if (interact) {
      Dud_count <- Dud_count + 1
      writeLines(as.character(Dud_count), Dud_fn)
    }
  } else {
    commfn <- file.path(comm_dn, paste0("node", indx, ".txt"))
    file.create(commfn)
    writeLines(as.character(B_new), commfn)
  }
  update_info <- list("mean_jaccards" = mean_jaccards, 
                      "consec_jaccards" = consec_jaccards,
                      "consec_sizes" = consec_sizes,
                      "found_cycle" = found_cycle,
                      "found_break" = found_break)
  
  # Checking overlap with previous sets
  if (interact && length(comms) > 0) {
    OL_check <- unlist(lapply(comms, function (C) jaccard(B_new, C)))
    if (sum(OL_check < 1 - OL_thres) > 0) {
      OL_count <- OL_count + 1
      writeLines(as.character(OL_count), OL_fn)
    }
  }
  
  # Noting which final.sets are filled; doing checks
  current_time <- proc.time()[3] - start_second
  if (interact && (current_time > time_limit || Dud_count > Dud_tol || OL_count > OL_tol)) {
    writeLines("TRUE", stop_fn)
  }
  
  return(list("StableComm" = B_new,
              "update_info" = update_info,
              "initial_set" = initial_set,
              "itCount" = itCount, "did_it_cycle" = did_it_cycle,
              "current_time" = current_time,
              "report" = "complete_extraction"))
}