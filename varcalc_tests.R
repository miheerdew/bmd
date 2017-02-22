source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")

m <- 20
my <- 10
rho <- 0.4
nsims <- 1000
ndata <- 1000
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- corsums_multi <- vars <- vars_multi <- rep(0, nsims * my)
simcorsums <- simvars <- rep(0, my)
timers <- timers_multi <- rep(0, nsims)

set.seed(12345)

for (sim in 1:nsims) {
  
  simpos <- (sim - 1) * my
  posindx <- (simpos + 1):(simpos + m)
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + matrix(rnorm(ndata * my, sd = sqrt(s2)),
                                         ncol = my)
  
  
  # Looped calculation

  timer <- get_nanotime()
  
  for (yindx in 1:my) {
    simcorsums[yindx] <- sum(cor(Yvars[ , yindx], Data))
    simvars[yindx] <- varcalc1(as.vector(Yvars[ , yindx]), Data)
  }
  corsums[posindx] <- simcorsums
  vars[posindx] <- simvars
  
  timers[sim] <- get_nanotime() - timer

  
  # Vectorized calculation
  
  timer <- get_nanotime()
  
  corsums_multi[posindx] <- as.vector(rowSums(cor(Yvars, Data)))
  vars_multi[posindx] <- varcalc1_multi(Yvars, Data)
  
  timers_multi[sim] <- get_nanotime() - timer
  
  
  if (sum((vars_multi[posindx] - vars[posindx])^2) > 1e-15 * my)
    break
  
}

if (!dir.exists('var_plots'))
  dir.create('var_plots')

png(file.path('var_plots', 'var_check_hom.png'))
hist(vars, main = '')
abline(v = popvar(rho), col = 'red', lwd = 2)
abline(v = mean(vars_multi), col = 'green', lwd = 2, lty = 2)
dev.off()


png(file.path('var_plots', 'z_check_hom.png'))
qqnorm(sqrt(ndata) * (corsums_multi - m * pxy(rho)) / sqrt(vars_multi))
abline(0, 1, col = 'red')
dev.off()

tmplot <- timers_multi / 1e9
tmplot <- tmplot[tmplot < quantile(tmplot, .99) & 
                   tmplot > quantile(tmplot, .01)]

png(file.path('var_plots', 'timer_hom.png'))
hist(tmplot, main = 'Timers for multiple Ys')
dev.off()

png(file.path('var_plots', 'timer_compare_hom.png'))
hist(log2(timers / timers_multi), main = '')
dev.off()

fivenum(timers / 1e9)
