source("sim_eQTL_network.R")
source("bmd.R")
total_expers <- readLines("sims-results/exper-names.txt")
run_expers <- c(10)

# This should consistent throughout the experiments
# (and match the same variable in sims/lfr/make_lfr_sims.R)
nreps <- 20

for (exper in run_expers) {
  
  exper_string <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", exper_string)
  
  # Loading parameters
  load(paste0(file.path("sims-results/sbm-par-lists", exper_string),
              ".RData"))
  
  for (p in 1:par_divs) {
    
    curr_dir_p <- file.path(root_dir, par_dirs[p])
    
    for (rep in 1:nreps) {
      
      cat("exper", exper, "p", p, "rep", rep, "\n")
      
      curr_dir_p_rep <- file.path(curr_dir_p, rep)
      load(file.path(curr_dir_p_rep, "sim.RData"))
      results <- bmd(sim$X, sim$Y, updateOutput = FALSE, OL_tol = 100, Dud_tol = 50)
      save(results, file = file.path(curr_dir_p_rep, "bmd.RData"))
        
      rm(sim, results)
      gc()
      
    }
    
  }
  
}