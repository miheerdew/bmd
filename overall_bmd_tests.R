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
res <- benchmark(C_F = bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE),
                 R_F_P = bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE, parallel = TRUE),
                 R_F = bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, calc_full_cor = TRUE),
                 C = bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE),
                 R_P = bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE),
                 R = bmd(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE),
                 C_P = bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE),
                 C_F_P = bmdC(sim$X, sim$Y, generalOutput = FALSE, verbose = FALSE, parallel = TRUE, calc_full_cor = TRUE),
                 replications = 2)
                
