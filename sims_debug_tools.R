source("sims_config.R")

#Size of all the X (or Y) variables
total_variables <- (nBg + nBM*nB)*sB
signal_nodes <- nBM*nB*sB

is_X <- function(u){
  u <= total_variables
}

bmd_index <- function(u){
  #Return the index of the bimodule u belongs to
  #Otherwise 0 if it's a noise node.

  if (u > total_variables){
    #So it must be a Y variable
    u <- u - total_variables
  }

  if( u <= signal_nodes ){
    #It is signal node. Return the index of the bimodule.
    return(1 + (u-1) %/% (nB*sB))
  }

  #Otherwise it must be noise.
  return(0)
}
