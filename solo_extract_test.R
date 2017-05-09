library(dplyr)
source("sim_eQTL_network.R")
source("auxiliary.R", local = TRUE)
source("bh_reject.R", local = TRUE)
source("pvals.R", local = TRUE)
source("initialize.R", local = TRUE)
source("extractNW.R", local = TRUE)
source("cbceNW.R",local=TRUE)

set.seed(1234567)



run <- function(data,initialBC, print_output=FALSE){
  updateOutput <<- TRUE
  alpha <<- 0.05
  updateMethod <<- 1
  Cpp <<- FALSE
  twoSided <<- FALSE
  calc_full_cor <<- TRUE
  # Set-up --------------------------
  start_second <<- proc.time()[3]
  X <<- scale(data$X); Y <<- scale(data$Y)
  X3 <<- X^3; X2 <<- X^2; X4ColSum <<- colSums(X^4)
  Y3 <<- Y^3; Y2 <<- Y^2; Y4ColSum <<- colSums(Y^4)

  dx <<- ncol(X)
  dy <<- ncol(Y)
  n  <<- nrow(X)

  Xindx <<- 1:dx
  Yindx <<- (dx + 1):(dx + dy)

  if(calc_full_cor){
    full_xy_cor <<- cor(X, Y)
  }
  extract(indx = 1, initialBC = initialBC, print_output = print_output)
}

experiment1 <- function(alphas, nsims, params){
  N <- length(alphas) * length(nsims)
  
  nsteps <- rep(NA,N)
  noise <- rep(list(NULL),N)
  missed <- rep(list(NULL),N)
  cycled <- rep(FALSE,N)
  alpha_seq <- rep(NA, N) 
  i <- 0
  
  for(n in 1:nsims){
   cat(sprintf("n=%d:\n alpha:",n))
   for(alpha in alphas){
      i <- i + 1
      cat(sprintf("%f,",alpha))
      data <- sim_eQTL_network(params)
      initialBC <- data$bms[[1]]
      fixedPnt <- try(run(data, initialBC,print_output=FALSE))
      if(inherits(fixedPnt, "try-error")){
        cat("Error occured")
        next
      }
      nsteps[i] <- fixedPnt$itCount
      noise[[i]] <- setdiff(fixedPnt$StableComm,initialBC)
      missed[[i]] <- setdiff(initialBC,fixedPnt$StableComm)
      cycled[i] <- fixedPnt$did_it_cycle
      alpha_seq[i] <- alpha
   }
  cat("\n")
  }
  data.frame(nsteps=nsteps,
            noise=I(noise),
            missed=I(missed),
            cycled=cycled,
            alphas = alpha_seq)
}

analyze <- function(result){
  result %>% group_by(alphas) %>% 
    summarise(
      mean_steps=mean(nsteps), 
      mean_noise=mean(sapply(noise,length)),
      mean_missed=mean(sapply(missed,length)),
      tot_cycled=sum(cycled)
    )
}

alphas <- c(0.00001, 0.00005 ,0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
nsims <- 500
#Large sample size, single bimodule.
param <- make_param_list(n = 1000, b = 1, cmin = 50, cmax = 100, bgmult = 2)
result <- experiment1(alphas,nsims,param)
save(result, file="solo_large.RData")
analysis <- analyze(result)
analysis

