source("sim_eQTL_network.R")

# For the experiments, do you want to shove the dec vec up to one?
shove_dec <- TRUE

# Global settings
par_divs <- 9
par_seq_dec <- 1:par_divs / (par_divs + 1)
par_seq  <- round(100 * (1:par_divs / (par_divs + 1)))
par_dirs <- as.character(par_seq)

# Give the names of your experiments: must be manually entered.
total_expers <- as.character(1:7)

if (!dir.exists("sims-results/sbm-par-lists"))
  dir.create("sims-results/sbm-par-lists", recursive = TRUE)

# Experiment 1 -----------------------------------------------------------------

main_text <- "Increase n"
par_list <- make_param_list()
pars <- c("n")
xlab <- "Sample Size"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(2000 * (par_seq_dec + min(par_seq_dec) * 
                                     as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment1.RData")

# Experiment 2 -----------------------------------------------------------------

main_text <- "Decrease min bimodule-half size"
par_list <- make_param_list()
pars <- c("cmin")
xlab <- "Min BM-half size"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- round(100 * (par_seq_dec + min(par_seq_dec) * 
                                    as.numeric(shove_dec)))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment2.RData")

# Experiment 3 -----------------------------------------------------------------

main_text <- "Increase amount of background"
par_list <- make_param_list()
pars <- c("bgmult")
xlab <- "#BG vars/#BM vars"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 10 * (par_seq_dec)

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment3.RData")

# Experiment 4 -----------------------------------------------------------------

main_text <- "Increase mean of beta parameters"
par_list <- make_param_list()
pars <- c("betamean")
xlab <- "Mean of beta params"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 10 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment4.RData")

# Experiment 5 -----------------------------------------------------------------

main_text <- "Increase eQTL probability"
par_list <- make_param_list()
pars <- c("p")
xlab <- "eQTL probability"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 1 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment5.RData")

# Experiment 6 -----------------------------------------------------------------

main_text <- "Increase intra-X correlation"
par_list <- make_param_list()
pars <- c("rho")
xlab <- "Intra-X corr"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 0.5 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment6.RData")

# Experiment 7 -----------------------------------------------------------------

main_text <- "Increase base noise scaling"
par_list <- make_param_list()
pars <- c("s2")
xlab <- "Noise scaling"
axis_par <- 1
par_settings <- matrix(0, 1, par_divs)
par_settings[1, ] <- 10 * (par_seq_dec + min(par_seq_dec) * as.numeric(shove_dec))

save(par_list,
     main_text,
     axis_par,
     pars,
     xlab,
     par_settings,
     par_seq,
     par_divs,
     par_dirs,
     file = "sims-results/sbm-par-lists/experiment7.RData")