library(lineprof)
source("bmd.R")
source("sim_eQTL_network.R")

par_list <- make_param_list(b=50)
sim <- sim_eQTL_network(par_list)
l <- lineprof(bmd(sim$X, sim$Y, calc_full_cor=TRUE))