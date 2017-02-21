nX <- 1
bX <- 20
rhos <- 0.4
Beta <- 1
s2 <- 1
nsims <- 100
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
pbr <- function (r) {Beta * (1 - r + bX * r)}
pxy <- function (r) {pbr(r) / sqrt(varY(r))}

# Function for single variance calculation
varY <- function (r) {bX * Beta * (1 - r + r * bX) + s2}
mu1111 <- function (r) {3}
mu1112 <- function (r) {3 * r}
mu1123 <- function (r) {2 * r^2 + r}
mu1122 <- function (r) {1 + 2 * r^2}
mu1234 <- function (r) {3 * r * (r - 1)^2 * (2 * r^2 + r) / d(r)}

# Functions for total variance calculation
pyy12  <- function (r) {((bX - 2) * mu1123(r) + (bX - 2) * (bX - 3) * mu1234(r) +
    4 * (bX - 2) * mu1123(r) + 2 * (mu1112(r) + mu1122(r)) + s2 * r) / varY(r)}
pyy11 <- function (r) {((bX - 1) * mu1122(r) + (bX - 1) * (bX - 2) * mu1123(r) +
    2 * (bX - 1) * mu1112(r) + mu1111(r) + s2) / varY(r)}
all4s <- function (r) {3 * sum(SigmaX)^2}
pyyyy <- function (r) {(all4s(r) + 
    6 * s2 * Beta^2 * bX * (1 - r + bX * r) + s2^2 * mu1111(r)) / varY(r)^(1.5)}
pyyjk <- function (r) {(2 * pbr(r)^2 * (2 * r + 1) / (r + 1) + 
                          r * Beta * pbr(r) * (bX + s2 / (Beta * pbr(r)) - 
                              2 * (1 - r + bX * r) / (1 + r))) / varY(r)}
pyyyy <- function (r) {3}
pyyjj <- function (r) {(2 * pbr(rhos)^2 + varY(rhos)) / varY(rhos)}
pyyyj <- function (r) {3 * varY(rhos) * pbr(rhos) / varY(rhos)^(3 / 2)}
pyjjk <- function (r) {pbr(r) * (2 * r + 1) / sqrt(varY(r))}
pyjjj <- function (r) {pbr(r) * 3 / sqrt(varY(r))}

# Variance calculations
popvar1234 <- function (r) {mu1234(r) + 0.25 * r^2 * 
    (mu1111(r) + 3 * mu1122(r)) - rho * 2 * mu1123(r)}

popvar1213 <- function (r) {mu1123(r) + 0.25 * r^2 * 
    (mu1111(r) + 3 * mu1122(r)) - rho * (mu1112(r) + mu1123(r))}

popvar <- function (r) {(pyyjj(r) + 0.25 * pxy(r)^2 * (pyyyy(r) + 2 * pyyjj(r) +
    mu1122(r)) - pxy(r) * (pyyyj(r) + pyjjj(r)) + 
    (bX - 1) * (pyyjk(r) + 0.25 * pxy(r)^2 * (pyyyy(r) + 2 * pyyjj(r) +
    mu1122(r)) - pxy(r) * (pyyyj(r) + pyjjk(r)))) * bX}

moments <- matrix(0, nsims, 5)
names(moments) <- c("mu1111", "mu1112", "mu1123", "mu1122", "mu1234")
moments <- as.data.frame(moments)
varYs <- rep(0, nsims)

crosscors <- matrix(0, nsims, 4)
crosscors <- as.data.frame(crosscors)
colnames(crosscors) <- c("ryy12", "ryy11", "ryyyy", "all4s")

truerhos <- matrix(0, nsims, 7)
truerhos <- as.data.frame(truerhos)
colnames(truerhos) <- c("cross", "pyyjk", "pyyyy", "pyyjj", "pyyyj", "pyjjk",
                        "pyjjj")

covars <- matrix(0, nsims, 4)
covars <- as.data.frame(covars)
colnames(covars) <- c("ysum", "r12", "r13", "r34")

for (sim in 1:nsims) {
  
  cat("sim", sim, "\n")
  
  Data <- mvrnorm(ndata, rep(0, bX), SigmaX)
  Y <- Data %*% rep(Beta, bX) + rnorm(ndata, sd = sqrt(s2))
  moments[sim, 'mu1111'] <- mean(Data[ , 1]^4)
  moments[sim, 'mu1112'] <- mean(Data[ , 1]^3 * Data[ , 2])
  moments[sim, 'mu1123'] <- mean(Data[ , 1]^2 * Data[ , 2] * Data[ , 3])
  moments[sim, 'mu1122'] <- mean(Data[ , 1]^2 * Data[ , 2]^2)
  moments[sim, 'mu1234'] <- mean(Data[ , 1] * Data[ , 2] * 
                                   Data[ , 3] * Data[ , 4])
  varYs[sim] <- var(Y)
  
  crosscors[sim, 'ryy12'] <- mean(Y^2 * Data[ , 1] * Data[ , 2]) /
    (var(Y) * sd(Data[ , 1] * sd(Data[ , 2])))
  crosscors[sim, 'ryy11'] <- mean(Y^2 * Data[ , 1]^2) / 
    (var(Y) * var(Data[ , 1]))
  crosscors[sim, 'ryyyy'] <- mean(Y^4) / sd(Y)^4
  crosscors[sim, 'all4s'] <- mean((Data %*% rep(1, bX))^4)
  
  sdY <- sd(Y); sdj <- sd(Data[ , 1]); sdk <- sd(Data[ , 2])
  truerhos[sim, 'cross'] <- mean(Y * Data[ , 1])
  truerhos[sim, 'pyyjk'] <- mean(Y^2 * Data[ , 1] * Data[ , 2]) / 
    (sdY^2 * sdj * sdk)
  truerhos[sim, 'pyyyy'] <- mean(Y^4) / (sdY^4)
  truerhos[sim, 'pyyjj'] <- mean(Y^2 * Data[ , 1]^2) / (sdY^2 * sdj^2)
  truerhos[sim, 'pyyyj'] <- mean(Y^3 * Data[ , 1]) / (sdY^3 * sdj)
  truerhos[sim, 'pyjjk'] <- mean(Y * Data[ , 1]^2 * Data[ , 2]) / 
    (sdY * sdj^2 * sdk)
  truerhos[sim, 'pyjjj'] <- mean(Y * Data[ , 1]^3 / (sdY * sdj^3))
  
  covars[sim, "ysum"] <- sum(cor(Y, Data))
  covars[sim, "r12"] <- cor(Data[ , 1], Data[ , 2])
  covars[sim, "r13"] <- cor(Data[ , 1], Data[ , 3])
  covars[sim, "r34"] <- cor(Data[ , 3], Data[ , 4])
  
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

png(file.path('moment_plots', 'moment_plot3.png'))
par(mfrow = c(3, 3))

hist(truerhos$cross)
abline(v = pbr(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$cros), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyyjk)
abline(v = pyyjk(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyyjk), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyyyy)
abline(v = pyyyy(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyyyy), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyyjj)
abline(v = pyyjj(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyyjj), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyyyj)
abline(v = pyyyj(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyyyj), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyjjk)
abline(v = pyjjk(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyjjk), col = 'green', lty = 2, lwd = 2)

hist(truerhos$pyjjj)
abline(v = pyjjj(rhos), col = 'red', lwd = 2)
abline(v = mean(truerhos$pyjjj), col = 'green', lty = 2, lwd = 2)

dev.off()

png(file.path('moment_plots', 'covars.png'))
par(mfrow = c(1, 3))

