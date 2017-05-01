source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVarMat.R")

m <- 20
rho <- 0.5
nsims <- 1000
ndata <- 1000
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- vars <- rep(0, nsims)

set.seed(12345)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  
  
  # variance calculation
  corsums[sim] <- sum(cor(Yvars, Data))
  vars[sim] <- varcalc1(Yvars, Data)
  
}
