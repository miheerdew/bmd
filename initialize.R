initialize1R <- function(u){
  
  if (u <= dx) {
    
    cors < ifelse(calc_full_cor, t(full_xy_cor[u,]), cor(Y, X[,u]))
    fischer_tranformed_cor <- atanh(cors) * sqrt(n - 3)
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
    successes <- Yindx[bh_reject(pvals, alpha)]
    
  } else {
    
    cors <- ifelse(calc_full_cor, full_xy_cor[ , u], cor(X, Y[ , u]))
    fischer_tranformed_cor <- atanh(cors) * sqrt(n - 3)
    pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
    successes <- Xindx[bh_reject(pvals, alpha)]

  }
  
  return(successes)
  
}


initialize1Cpp <- function(u) {
  
  if (u > dx) {
    
    #Test X
    u <- u - dx
    cors <- ifelse(calc_full_cor, full_xy_cor[ , u], cor(X, Y[ , u]))
    return(initializeC(n, cors, alpha, conserv))
    
  } else {
    
    #Test Y
    cors <- ifelse(calc_full_cor, t(full_xy_cor[u,]), cor(Y, X[,u]))
    return(initializeC(n, cors, alpha, conserv) + dx)
    
  }
  
}