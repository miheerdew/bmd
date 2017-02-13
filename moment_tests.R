nX <- 1
bX <- 20
rhos <- 0.5
Beta <- 1
s2 <- 1
nsims <- 100
ndata <- 10000

library(MASS)
library(Matrix)

# Set intra-correlations of X's
rho_blocksX <- lapply(rhos, function (R) matrix(R, bX, bX) + 
                        diag(rep(1 - R, bX)))
mX <- nX * bX

#------Generating Block Covariance Matrices-----#
SigmaX <- as.matrix(bdiag(rho_blocksX))

# Useful functions
d <- function (r) {1 - 3 * r^2 + 2 * r^3}
varY <- function (r) {bX * Beta * (1 - r + r * bX) + s2}
mu1111 <- function (r) {3}
mu1112 <- function (r) {3 * r}
mu1123 <- function (r) {2 * r^2 + r}
mu1122 <- function (r) {1 + 2 * r^2}
mu1234 <- function (r) {3 * r * (r - 1)^2 * (2 * r^2 + r) / d(r)}
ryy12  <- function (r) {(bX - 2)^2 * mu1234(r) + 4 * (bX - 2) * mu1123(r) +
    2 * (mu1122(r) + mu1112(r)) + s2 * r}
ryy11 <- function (r) {(bX - 2)^2 * mu1123(r) + 4 * (bX - 2) * mu1112(r) +
    4 * mu1111(r) + s2}

moments <- matrix(0, nsims, 5)
names(moments) <- c("mu1111", "mu1112", "mu1123", "mu1122", "mu1234")
moments <- as.data.frame(moments)
varYs <- rep(0, nsims)

crosscors <- matrix(0, nsims, 2)
names(crosscors) <- c("ryy12", "ryy11")
crosscors <- as.data.frame(crosscors)

for (sim in 1:nsims) {
  
  cat("sim", sim, "\n")
  
  Data <- mvrnorm(ndata, rep(0, bX), SigmaX)
  Y <- Data %*% rep(1, bX) + rnorm(ndata, sd = sqrt(s2))
  moments[sim, 'mu1111'] <- mean(Data[ , 1]^4)
  moments[sim, 'mu1112'] <- mean(Data[ , 1]^3 * Data[ , 2])
  moments[sim, 'mu1123'] <- mean(Data[ , 1]^2 * Data[ , 2] * Data[ , 3])
  moments[sim, 'mu1122'] <- mean(Data[ , 1]^2 * Data[ , 2]^2)
  moments[sim, 'mu1234'] <- mean(Data[ , 1] * Data[ , 2] * Data[ , 3] * Data[ , 4])
  varYs[sim] <- var(Y)
  
  crosscors[sim, 'ryy12'] <- mean(Y^2 * Data[ , 1] * Data[ , 2])
  crosscors[sim, 'ryy11'] <- mean(Y^2 * Data[ , 1]^2)
  
}

if (!dir.exists('moment_plots'))
  dir.create('moment_plots')

png(file.path('moment_plots', 'moment_plot1.png'))
par(mfrow = c(2, 3))

hist(moments$mu1111)
abline(v = mu1111(rhos), col = 'red')

hist(moments$mu1112)
abline(v = mu1112(rhos), col = 'red')

hist(moments$mu1122)
abline(v = mu1122(rhos), col = 'red')

hist(moments$mu1123)
abline(v = mu1123(rhos), col = 'red')

hist(moments$mu1234)
abline(v = mu1234(rhos), col = 'red')

hist(varYs)
abline(v = varY(rhos), col = 'red')

dev.off()

png(file.path('moment_plots', 'moment_plot2.png'))
par(mfrow = c(1, 2))

hist(crosscors$ryy12)
abline(v = ryy12(rhos), col = 'red')
abline(v = mean(crosscors$ryy12), col = 'green')

hist(crosscors$ryy11)
abline(v = ryy11(rhos), col = 'red')
abline(v = mean(crosscors$ryy11), col = 'green')

dev.off()