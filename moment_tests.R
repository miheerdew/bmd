nX <- 1
bX <- 20
rhos <- 0.5
Beta <- 1
s2 <- 1

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
varY <- function (r) {bX * Beta * (1 - rho + rho * bX) + s2}
mu1111 <- function (r) {3}
mu1112 <- function (r) {3 * rho}
mu1123 <- function (r) {2 * rho^2 + rho}
mu1122 <- function (r) {1 + 2 * rho^2}
mu1234 <- function (r) {3 * rho * (rho - 1)^2 * (2 * rho^2 + rho) / d(p)}

S4 <- diag(4)
S4[col(S4) != row(S4)] <- rho
S3 <- S4[1:3, 1:3]
S3[col(S3) != row(S3)] <- rho
one3 <- rep(1, 3)
t(one3) %*% solve(S3) %*% one3
3 * (rho - 1)^2 / d(rho)

e <- function (r) { 1 - 3 * r^2 * (r - 1)^2 / d(r) }

rhos <- seq(0.01, 0.99, b = 0.01)
plot(rhos, sapply(rhos, e))


library(MASS)
S2 <- S3[1:2, 1:2]
Data2 <- mvrnorm(10000, rep(0, 2), S2)
mean(Data2[ , 1]^2 * Data2[ , 2]^2)
1 + 2 * rho^2


Data3 <- mvrnorm(10000, rep(0, 3), S3)
mean(Data3[ , 1]^2 * Data3[ , 2] * Data3[ , 3])

2 * rho^2 + rho

Data4 <- mvrnorm(10000, rep(0, 4), S4)
mean(apply(Data4, 1, prod))
3* rho * (rho - 1)^2 * (2 * rho^2 + rho) / d(rho)