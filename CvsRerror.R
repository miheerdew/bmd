load("CvsR/CvsR.RData")

X3 <- X^3; X2 <- X^2; X4ColSum <- colSums(X^4)
Y3 <- Y^3; Y2 <- Y^2; Y4ColSum <- colSums(Y^4)

dx <- ncol(X)
dy <- ncol(Y)
n  <- nrow(X)

Xindx <- 1:dx
Yindx <- (dx + 1):(dx + dy)
calc_full_cor <- FALSE
alpha <- 0.05

library(bmdCpp)
source("bh_reject.R", local = TRUE)
source("pvals.R", local = TRUE)


Cpvals <- c(pvalsCpp(B_oldy), pvalsCpp(B_oldx))
Rpvals <- c(pvalsR(B_oldy), pvalsR(B_oldx))


CB_new <- bh_rejectC(Cpvals, alpha, conserv = TRUE)
RB_new <- bh_rejectR(Rpvals, alpha)

identical(CB_new, RB_new)
# False
