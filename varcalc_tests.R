library(MASS)

source("varcalcs.R")
source("moment_calcs.R")

m <- 20
rho <- 0.4
nsims <- 1000
ndata <- 1000
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- vars <- rep(0, nsims)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnorm(ndata, rep(0, m), SigmaX_hom)
  Y <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  corsums[sim] <- sum(cor(Y, Data))
  vars[sim] <- varcalc1(Data, as.vector(Y))
  
}