library(MASS)
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")

m <- 20
rho <- 0.4
nsims <- 1000
ndata <- 1000
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- vars <- timers <- rep(0, nsims)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnorm(ndata, rep(0, m), SigmaX_hom)
  Y <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  corsums[sim] <- sum(cor(Y, Data))
  
  timer <- get_nanotime()
  vars[sim] <- varcalc1(Data, as.vector(Y))
  timers[sim] <- get_nanotime() - timer
  
}

if (!dir.exists('var_plots'))
  dir.create('var_plots')

png(file.path('var_plots', 'var_check_hom.png'))
hist(vars, main = '')
abline(v = popvar(rho), col = 'red', lwd = 2)
abline(v = mean(vars), col = 'green', lwd = 2, lty = 2)
dev.off()


png(file.path('var_plots', 'z_check_hom.png'))
qqnorm(sqrt(ndata) * (corsums - m * pxy(rho))/sqrt(vars))
abline(0, 1, col = 'red')
dev.off()

png(file.path('var_plots', 'timer_hom.png'))
hist(timers / 1e9, main = '')
dev.off()

fivenum(timers / 1e9)
