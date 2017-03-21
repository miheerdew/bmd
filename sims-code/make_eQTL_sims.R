source("sim_eQTL_network.R")

total_expers <- readLines("sims-results/exper-names.txt")

run_expers <- c(3, 8:10)

redrawSeeds <- FALSE

# This should consistent throughout the experiments
nreps <- 20

for (exper in run_expers) {
  
  set.seed(12345)
  
  exper_string <- paste0("experiment", total_expers[exper])
  
  # Finding the folder
  root_dir <- file.path("sims-results", exper_string)
  if (!dir.exists(root_dir)) {dir.create(root_dir)}
  
  # Loading the parameters (change them with sims/sbm_sims/make_par_lists.R)
  load(paste0("sims-results/sbm-par-lists/", exper_string, ".RData"))
  
  # Getting par settings
  par_divs <- ncol(par_settings)
  par_points <- 1:par_divs / (par_divs + 1)
  
  # Getting rep folder names
  rep_dirs <- as.character(1:nreps)
  
  
  for (p in 1:par_divs) {
    
    cat("Making sim", exper_string, "par num", p, "- ")
    
    
    # Setting directory
    curr_dir_p <- file.path(root_dir, par_dirs[p])
    if (!dir.exists(curr_dir_p)) {dir.create(curr_dir_p)}
    
    # Getting par_list_p
    par_list_p <- par_list
    for (j in 1:length(pars)) {
      par <- pars[j]
      par_indx <- which(names(par_list) == par)
      par_list_p[[par_indx]] <- par_settings[j, p]
    }
    
    
    
    for (rep in 1:nreps) {
      
      if (rep == 1)
        cat("rep: ")
      cat(rep)
      if (rep == nreps)
        cat("\n")
      
      # Setting directory
      curr_dir_p_rep <- file.path(curr_dir_p, rep_dirs[rep])
      if (!dir.exists(curr_dir_p_rep)) {dir.create(curr_dir_p_rep)}
      sim_fn <- file.path(curr_dir_p_rep, "sim.RData")
      
        
      # Draw random seed and save
      seedfn <- file.path(curr_dir_p_rep, "sim_seed.txt")
      if (redrawSeeds || !file.exists(seedfn)) {
        seed_draw <- sample(1e6, 1)
        writeLines(as.character(seed_draw), seedfn)
        set.seed(seed_draw)
      } else {
        set.seed(as.numeric(readLines(file.path(curr_dir_p_rep, 
                                                "sim_seed.txt")
        )))
      }
      
      sim <- sim_eQTL_network(par_list_p)
      save(sim, file = sim_fn)
      rm(sim)
      gc()
        
    } 
    
  }
  
}