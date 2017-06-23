source("cbceNW.R")
library(foreach)
library(doParallel)
library(doRNG)


cores=detectCores()-1
cl <- makeCluster(cores, output="type1error.txt")
clusterEvalQ(cl, source("cbceNW.R"))
registerDoParallel(cl)

ns <- 1000#seq(0,1000,by=50)
alphas <- 0.05#seq(0,0.5, 0.001)
m <- 500
nsims <- 100

set.seed(12345)

generate_data <- function(dx, dy, m){
 X <- matrix(rnorm(dx*m), ncol=dx)
 Y <- matrix(rnorm(dy*m), ncol=dy)
 list(X=X, Y=Y)
}

type1error <- function(dx, dy, m, alpha){
  #Compute type 1 error by the above params by simulation
  r <- foreach(s = 1:nsims,
                 .combine = c) %do% {
            data <- generate_data(dx, dy, m)
            res <-
              cbceNW(
                data$X,
                data$Y,
                alpha = alpha,
                exhaustive = TRUE,
                calc_full_cor = TRUE,
                Cpp = FALSE
              )
            length(res$finalIndxs) > 0
          }
  mean(r)
}

N <- length(alphas) * length(ns)
res <- data.frame(alpha=rep(NA, N), n=rep(NA, N), type1error=rep(NA,N))
i <- 1

for(alpha in alphas){
  for(n in ns){
   res[i,]  <- list(alpha, n, type1error(n, n, m, alpha))
   i <- i + 1
  }
}

stopCluster(cl)