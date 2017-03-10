source("sims_config.R")
source("sims_debug_tools.R")
source("bmd.R")


n <- 3000

tracer <- function(){
  if(deparse(sys.call(-5)) == "bh_reject(pvals, alpha)"){
    f <- sys.frame(-5)
    plot_pvals_by_bmd(f$pvals, f$alpha/f$mults, transformed=TRUE,cap=100)
    x <- readline(prompt="Press [enter] to continue;")
    if (x != "") {
      stop("interrupted")
    }
  }
}

################
# Start execution
setBreakpoint("bmd.R#130", tracer = tracer)
load(dataset_fname(n))
BMDresults <- bmd(X, Y, tag = n,
                  saveDir = file.path(saveDir, "BMD_saves"),
                  updateMethod = 5, initializeMethod = 3,
                  Dud_tol = 10, OL_tol = 10, time_limit = 1800,
                  bmd_index = bmd_index, calc_full_cor=TRUE)