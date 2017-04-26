source("cbceNW.R")
source("cbce.R")
load("sims-results/experiment1/90/1/sim.RData")
library(rbenchmark)

res <- benchmark(NW = cbceNW(sim$X, sim$Y, calc_full_cor = TRUE, verbose = FALSE),
                 W  = cbce(sim$X, sim$Y, calc_full_cor = TRUE, verbose = FALSE),
                 replications = 10)
res

