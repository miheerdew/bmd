# Compares two clusterings of bimodules.
# Clusterings must have the following format:
#   - List object of two items:
#       1. List object for communities:
#           a. List object for X variable sets
#           b. List object for Y variable sets
#               * must be same length as X var sets
#       2. Vector object contating background indices
#
# Computes the following metrics:
#   - Best match distance between communities
#       * Returns 0 if no communities
#   - Jaccard distance between backgrounds
#-------------------------------------------------------------------------------
symdiff <- function (s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}

jaccard <- function (s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}

best_match_bimodule <- function (C1, C2, bg1, bg2) {
  
  n <- length(C1)
  m <- length(C2)
  J <- J1 <- J2 <- matrix(0, n, m)
  
  # Computing jaccard distance between background declarations
  
  if (length(bg1) + length(bg2) == 0) {
    
    BJ <- 0
    
  } else {
    
    BJ <- jaccard(bg1, bg2)
    
  }
  
  # Computing best match metric for communities
  
  if (n * m == 0) {
    
    BM <- 0
    
  } else {
    
    for (i in 1:n) {
      
      C1i <- C1[[i]]
      
      for (j in 1:m) {
        
        C2j <- C2[[j]]
        
        J[i, j] <- jaccard(C1i, C2j)
        J1[i, j] <- 1 - length(intersect(C1i, C2j)) / length(C1i)
        J2[i, j] <- 1 - length(intersect(C1i, C2j)) / length(C2j)
        
      }
      
    }
    
    BM <- (sum(apply(J, 1, min)) + sum(apply(J, 2, min))) / (2 * (n + m))
    BM1 <- (sum(apply(J1, 1, min)) + sum(apply(J1, 2, min))) / (2 * (n + m))
    BM2 <- (sum(apply(J2, 1, min)) + sum(apply(J2, 2, min))) / (2 * (n + m))
    

  }
  
  return(c("BestMatch" = 1 - BM, 
           "BestMatch1" = 1 - BM1,
           "BestMatch2" = 1 - BM2, 
           "BackgroundMatch" = 1 - BJ))
  
}
    
  
  
  
  
  
  