run_expers <- sapply(commandArgs(TRUE), as.numeric)

source("sim_eQTL_network.R")
source("bmd.R")
source("bmdC.R")
source("bmd_cpp.R")
source("run_brim.R")
source("ircc.R")
total_expers <- readLines("sims-results/exper-names.txt")

runBMDcpp <- TRUE
runBMD2 <- TRUE
runBMD <- FALSE
runBRIM <- FALSE
runkmeans <- FALSE

# This should consistent throughout the experiments
# (and match the same variable in sims/lfr/make_lfr_sims.R)
nreps <- 10

for (exper in run_expers) {
  
  exper_string <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", exper_string)
  
  # Loading parameters
  load(paste0(file.path("sims-results/sbm-par-lists", exper_string), ".RData"))
  
  for (p in 1:par_divs) {
    
    curr_dir_p <- file.path(root_dir, par_dirs[p])
    
    # Setting alpha
    alpha <- ifelse(palpha, par_settings[1, p], 0.05)
    
    for (rep in 1:nreps) {
      
      cat("exper", exper, "p", p, "rep", rep, "\n")
    
      curr_dir_p_rep <- file.path(curr_dir_p, rep)
      load(file.path(curr_dir_p_rep, "sim.RData"))
      
      # Running BMDcpp
      if (runBMDcpp) {
        timer <- proc.time()[3]
        results <- bmdC(sim$X, sim$Y, alpha = alpha,
                           updateOutput = FALSE, OL_tol = 100, Dud_tol = 50,
                           calc_full_cor = TRUE, updateMethod = 5)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd_cpp.RData"))
      }
      
      # Running BMD2
      if (runBMD2) {
        timer <- proc.time()[3]
        results <- bmd(sim$X, sim$Y, alpha = alpha,
                       updateOutput = FALSE, OL_tol = 100, Dud_tol = 50,
                       calc_full_cor = TRUE, updateMethod = 5)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd2.RData"))
      }
      
      # Running BMD
      if (runBMD) {
        timer <- proc.time()[3]
        results <- bmd(sim$X, sim$Y, alpha = alpha,
                       updateOutput = FALSE, OL_tol = 100, Dud_tol = 50,
                       calc_full_cor = TRUE, updateMethod = 7)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "bmd.RData"))
      }
      
      # Running BRIM
      if (runBRIM) {
        timer <- proc.time()[3]
        results <- run_brim(sim$X, sim$Y, alpha = alpha)
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "brim.RData"))
      }
      
      # Running IRCC-kmeans
      if (runkmeans) {
        nbmds <- ifelse(par_list$bgmult > 0, par_list$b + 1, par_list$b)
        timer <- proc.time()[3]
        results <- ircc(sim$X, sim$Y, nbmds = nbmds, method = "kmeans")
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir_p_rep, "kmeans.RData"))
      }
        
      rm(sim, results)
      gc()
      
    }
    
  }
  
}