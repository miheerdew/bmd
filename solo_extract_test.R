source("sim_eQTL_network.R")
simple_net_pars <- make_param_list(n = 200, b = 1, cmin = 100, cmax = 100, bgmult = 1)
simple_net <- sim_eQTL_network(simple_net_pars)

# Set-up --------------------------
calc_full_cor <- TRUE
twoSided <- FALSE
Cpp <- FALSE
updateMethod <- 1
alpha <- 0.05
start_second <- proc.time()[3]
X <- scale(simple_net$X); Y <- scale(simple_net$Y)
X3 <- X^3; X2 <- X^2; X4ColSum <- colSums(X^4)
Y3 <- Y^3; Y2 <- Y^2; Y4ColSum <- colSums(Y^4)

dx <- ncol(X)
dy <- ncol(Y)
n  <- nrow(X)

Xindx <- 1:dx
Yindx <- (dx + 1):(dx + dy)

if(calc_full_cor){
  full_xy_cor <- cor(X, Y)
}

# To assess whether or not the bi-module is a fixed point, you can run the "extract" function:
source("auxiliary.R", local = TRUE)
source("bh_reject.R", local = TRUE)
source("pvals.R", local = TRUE)
source("initialize.R", local = TRUE)
source("extractNW.R", local = TRUE)

initialBC <- c(1:100, 201:300)
fxdPnt <- extract(indx = 1, initialBC = initialBC, print_output = FALSE)
fxdPnt$StableComm
