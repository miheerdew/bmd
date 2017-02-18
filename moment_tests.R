nX <- 1
bX <- 20
rhos <- 0.4
Beta <- 1
s2 <- 1
nsims <- 1000
ndata <- 100000

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
pyy12  <- function (r) {((bX - 2) * mu1123(r) + (bX - 2) * (bX - 3) * mu1234(r) +
    4 * (bX - 2) * mu1123(r) + 2 * (mu1112(r) + mu1122(r)) + s2 * r) / varY(r)}
pyy11 <- function (r) {((bX - 1) * mu1122(r) + (bX - 1) * (bX - 2) * mu1123(r) +
    2 * (bX - 1) * mu1112(r) + mu1111(r) + s2) / varY(r)}
all4s <- function (r) {3 * sum(SigmaX)^2}
pyyyy <- function (r) {(all4s(r) + 
    6 * s2 * Beta^2 * bX * (1 - r + bX * r) + s2^2 * mu1111(r)) / varY(r)^(1.5)}

moments <- matrix(0, nsims, 5)
names(moments) <- c("mu1111", "mu1112", "mu1123", "mu1122", "mu1234")
moments <- as.data.frame(moments)
varYs <- rep(0, nsims)

crosscors <- matrix(0, nsims, 4)
names(crosscors) <- c("ryy12", "ryy11", "ryyyy", "all4s")
crosscors <- as.data.frame(crosscors)

for (sim in 1:nsims) {
  
  cat("sim", sim, "\n")
  
  Data <- mvrnorm(ndata, rep(0, bX), SigmaX)
  Y <- Data %*% rep(Beta, bX) + rnorm(ndata, sd = sqrt(s2))
  moments[sim, 'mu1111'] <- mean(Data[ , 1]^4)
  moments[sim, 'mu1112'] <- mean(Data[ , 1]^3 * Data[ , 2])
  moments[sim, 'mu1123'] <- mean(Data[ , 1]^2 * Data[ , 2] * Data[ , 3])
  moments[sim, 'mu1122'] <- mean(Data[ , 1]^2 * Data[ , 2]^2)
  moments[sim, 'mu1234'] <- mean(Data[ , 1] * Data[ , 2] * Data[ , 3] * Data[ , 4])
  varYs[sim] <- var(Y)
  
  crosscors[sim, 'ryy12'] <- mean(Y^2 * Data[ , 1] * Data[ , 2]) /
    (var(Y) * sd(Data[ , 1] * sd(Data[ , 2])))
  crosscors[sim, 'ryy11'] <- mean(Y^2 * Data[ , 1]^2) / (var(Y) * var(Data[ , 1]))
  crosscors[sim, 'ryyyy'] <- mean(Y^4) / sd(Y)^4
  crosscors[sim, 'all4s'] <- mean((Data %*% rep(1, bX))^4)
  
}

if (!dir.exists('moment_plots'))
  dir.create('moment_plots')

png(file.path('moment_plots', 'moment_plot1.png'))
par(mfrow = c(2, 3))

hist(moments$mu1111, main = '')
abline(v = mu1111(rhos), col = 'red', lwd = 2)
abline(v = mean(moments$mu1111), col = 'green', lty = 2, lwd = 2)

hist(moments$mu1112, main = '')
abline(v = mu1112(rhos), col = 'red', lwd = 2)
abline(v = mean(moments$mu1112), col = 'green', lty = 2, lwd = 2)

hist(moments$mu1122, main = '')
abline(v = mu1122(rhos), col = 'red', lwd = 2)
abline(v = mean(moments$mu1122), col = 'green', lty = 2, lwd = 2)

hist(moments$mu1123, main = '')
abline(v = mu1123(rhos), col = 'red', lwd = 2)
abline(v = mean(moments$mu1123), col = 'green', lty = 2, lwd = 2)

hist(moments$mu1234, main = '')
abline(v = mu1234(rhos), col = 'red', lwd = 2)
abline(v = mean(moments$mu1234), col = 'green', lty = 2, lwd = 2)

hist(varYs, main = '')
abline(v = varY(rhos), col = 'red', lwd = 2)
abline(v = mean(varYs), col = 'green', lty = 2, lwd = 2)

dev.off()

png(file.path('moment_plots', 'moment_plot2.png'))
par(mfrow = c(2, 3))

hist(crosscors$ryy12)
abline(v = pyy12(rhos), col = 'red', lwd = 2)
abline(v = mean(crosscors$ryy12), col = 'green', lty = 2, lwd = 2)

hist(crosscors$ryy11)
abline(v = pyy11(rhos), col = 'red', lwd = 2)
abline(v = mean(crosscors$ryy11), col = 'green', lty = 2, lwd = 2)

hist(crosscors$ryyyy)
abline(v = pyyyy(rhos), col = 'red', lwd = 2)
abline(v = mean(crosscors$ryyyy), col = 'green', lty = 2, lwd = 2)

hist(crosscors$all4s)
abline(v = all4s(rhos), col = 'red', lwd = 2)
abline(v = mean(crosscors$all4s), col = 'green', lty = 2, lwd = 2)

dev.off()