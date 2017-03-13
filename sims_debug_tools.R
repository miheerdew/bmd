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

bmd_proportion <- function(nodes){
    table(sapply(nodes, bmd_index))
}

#The indices corresponding to the bimodules. Accessible as list
bmd_indices <- c(lapply(1:nBM, function(b) 1:(nB*sB) + (b-1)*nB*sB))
noise_indices <- 1:(nBg*sB) + signal_nodes

plot_pvals_by_bmd <- function(pvals, alpha, ...){
  old.par <- par(mfrow=c(2,2), mar=c(2,2,2,2))
  for (i in 1:nBM){
    plot_pvals(pvals[bmd_indices[[i]]], alpha,
              bmd_name=paste("BMD",i), ...)
  }
  par(old.par)
}

plot_pvals <- function(pvals, alpha, bmd_name,
                       transformed=FALSE, cap=Inf){
  if(transformed){
    transform <- function(x) { -log10(x) }
    legend_loc <- "topright"
    ylab <- "-log(p(i))"
  } else {
    transform <- identity
    legend_loc <- "topleft"
    ylab <- "p(i)"
  }

  n <- length(pvals)
  pvals <- pmin(transform(sort(pvals)),cap)
  thresh_line  <- transform(c(1:n) * alpha / n)
  plot(pvals, xlab="i (p-val rank)", ylab=ylab, pch=1)
  lines(1:n, thresh_line, col = "purple")
  title(bmd_name)
}

plot_pval_simple <- function(pvals, alpha, group_index){
  n <- length(pvals)
  groups = sapply(order(pvals),group_index)
  group_nums = sort(unique(groups))
  plot(sort(pvals), col = 1 + groups)
  lines(1:n, c(1:n) * alpha / n, col = "purple")
  legend("topleft",
         legend = c(paste0("group ", group_nums), "BH line"),
         col = c(1 + group_nums, "purple"),
         lty = c(rep(NA, length(group_nums)), 1),
         pch = c(rep(1, length(group_nums)), NA)
        )
}

plot_pval_transformed <- function(pvals, alpha, group_index){
  n <- length(pvals)
  groups = sapply(order(pvals),group_index)
  group_nums = sort(unique(groups))
  plot(-log10(sort(pvals)), col = 1 + groups)
  lines(1:n, -log10(c(1:n) * alpha / n), col = "purple")
  legend("topright",
         legend = c(paste0("group ", group_nums), "BH line"),
         col = c(1 + group_nums, "purple"),
         lty = c(rep(NA, length(group_nums)), 1),
         pch = c(rep(1, length(group_nums)), NA)
        )
}
