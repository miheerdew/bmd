initialize1R <- function(u) {
  
  if (u <= dx) {
    
    # Test X
    cors < if (calc_full_cor) { 
      t(full_xy_cor[u,])
    } else {
      cor(Y, X[,u])
    }
    fischer_tranformed_cor <- atanh(cors) * sqrt(n - 3)
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
    successes <- Yindx[bh_reject(pvals, alpha)]
    
  } else {
    
    # Test Y
    cors <- if (calc_full_cor) {
      full_xy_cor[ , u]
    } else {
      cor(X, Y[ , u])
    }
    fischer_tranformed_cor <- atanh(cors) * sqrt(n - 3)
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
    successes <- Xindx[bh_reject(pvals, alpha)]

  }
  
  return(successes)
  
}


initialize1Cpp <- function(u) {
  
  if (u > dx) {
    
    # Test X
    u <- u - dx
    if (calc_full_cor) {
      cors <-  full_xy_cor[ , u]
    } else {
      cors <- cor(X, Y[ , u])
    }
    return(initializeC(n, cors, alpha, conserv = TRUE))
    
  } else {
    
    # Test Y
    cors <- if (calc_full_cor) {
      cors <- t(full_xy_cor[u, ]) 
    } else {
      cors <-  cor(Y, X[ , u])
    }
    return(initializeC(n, cors, alpha, conserv = TRUE) + dx)
    
  }
  
}