load("sims-results/experiment1/90/1/sim.RData")
source("bmdC.R")
source("bmd.R")
library(rbenchmark)

# Running non-parallel methods
CC <- bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE)
RC <- bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE)
C <- bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE)
R <- bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE)

# Formatting for equality check
CC_mod <- lapply(CC$extract_res, function (R) c(R$StableComm, R$initial_set))
RC_mod <- lapply(RC$extract_res, function (R) c(R$StableComm, R$initial_set))
C_mod <- lapply(C$extract_res, function (R) c(R$StableComm, R$initial_set))
R_mod <- lapply(R$extract_res, function (R) c(R$StableComm, R$initial_set))

# Equality check
identical(CC_mod, RC_mod); identical(RC_mod, C_mod); identical(R_mod, C_mod)

# Speed check; this will take about 10 minutes
res <- benchmark(bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE),
                 bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE, parallel = TRUE),
                 bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE),
                 bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE),
                 bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE),
                 bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE),
                 bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE),
                 bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE, calc_full_cor = TRUE),
                 replications = 5)
                
