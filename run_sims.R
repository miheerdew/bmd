source("sims_config.R")

RUN <- list(sims    =FALSE,
            methods =FALSE,
            plots   =TRUE)

if(RUN$sims){
  unlink(file.path(saveDir,"datasets"), recursive = TRUE)
  source("sims.R")
}

if(RUN$methods){
  unlink(file.path(saveDir,"results"), recursive = TRUE)
  source("sims_run_methods.R")
}

if(RUN$plots){
  unlink(file.path(saveDir,"plots"), recursive = TRUE)
  source("simplots.R")
}
