source("bmd_cpp.R")
source("sim_eQTL_network.R")

set.seed(123456)
par_list <- make_param_list(b=100)
sim <- sim_eQTL_network(par_list)
prof_file <- 'bmdprofiling.dat'
Rprof(prof_file)
bmd_cpp(sim$X, sim$Y, calc_full_cor=TRUE)
Rprof()
summaryRprof(prof_file)

