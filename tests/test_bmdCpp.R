library(testthat)
library(microbenchmark)
library(bmdCpp)

source("../varcalcs.R", chdir = T)
source("../sim_eQTL_network.R", chdir = T)

context("BMD helpers")

sim <- sim_eQTL_network(make_param_list(b = 20))
X <- scale(sim$X)
Y <- scale(sim$Y)
dy <- ncol(Y)
dx <- ncol(X)
n <- nrow(X)
alpha = 0.05


which_side <- function (A, test_x) {
  if (test_x) {
    return(list(second = Y[, A, drop=FALSE],
                first = X))
  } else {
    return(list(second = X[, A, drop=FALSE],
                first = Y))
  }
}

pvalsR <- function(firstMat, secondMat) {
  vars <- varcalc1_multi(firstMat, secondMat)
  zstats <- sqrt(n) * rowSums(cor(firstMat, secondMat)) / sqrt(vars)
  as.vector(pnorm(zstats, lower.tail = FALSE))
}

compare_pvalues <- function(A, test_x, benchmark = FALSE, ...) {
  mat <- which_side(A, test_x)
  f <- ifelse(benchmark, microbenchmark, expect_equivalent)
  f(
    pvalsR(mat$first, mat$second),
    pvalsC(
      mat$first,
      mat$second,
      colSums(mat$first ^ 4),
      mat$first ^ 2,
      mat$first ^ 3,
      cor(mat$first, mat$second)
    )
  )
}

test_that("Computes pvalues correctly", {
  compare_pvalues(1:3, FALSE)
  compare_pvalues(1:3, TRUE)
})

bh_rejectR <- function (pvals, alpha, conserv) {
  m <- length(pvals)

  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals)
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals)
  }

  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

compare_bhreject <-
  function(pval, alpha, conserv, benchmark = FALSE, ...) {
    f <- ifelse(benchmark, microbenchmark, expect_equivalent)
    f(bh_rejectR(pval, alpha, conserv),
      bh_rejectC(pval, alpha, conserv),
      ...)
  }

test_that("BH reject works correctly", {
  pval <- runif(10000)
  compare_bhreject(pval, 0.05, FALSE)

  pval <- c(runif(10000), 0, 1, 0)
  compare_bhreject(pval, 0.05, TRUE)

  pval <- c(runif(10000), 0, 1, 0) / 10
  compare_bhreject(pval, 0.05, TRUE)
})

initializeR <- function(cor_to_u, alpha, conserv) {
  fischer_tranformed_cor <- atanh(cor_to_u) * sqrt(n - 3)
  pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
  return(bh_rejectR(pvals, alpha, conserv))
}

compare_inititalization <- function(u,
                                   test_x,
                                   conserv = TRUE,
                                   benchmark = FALSE,
                                   ...) {
    if (test_x) {
      cor_to_u <- cor(X, Y[, u])
    } else {
      cor_to_u <- cor(Y, X[, u])
    }

    f <- ifelse(benchmark, microbenchmark, expect_equivalent)

    f(
      initializeR(cor_to_u, alpha, conserv),
      initializeC(n, cor_to_u, alpha, conserv),
      ...
    )
  }

test_that("initialize works correctly", {
  compare_inititalization(23, TRUE, TRUE)

  compare_inititalization(52, FALSE, TRUE)

  compare_inititalization(83, FALSE, FALSE)

  compare_inititalization(122, TRUE, FALSE)
})

compare_updates <-
  function(A, test_x, conserv, benchmark = FALSE, ...) {
    mat <- which_side(A, test_x)
    f <- ifelse(benchmark, microbenchmark, expect_equivalent)
    f4sum <-colSums(mat$first ^ 4)
    f2 <- mat$first ^ 2
    f3 <- mat$first ^ 3
    cross <- cor(mat$first, mat$second)

    f(
      bh_rejectR(pvalsR(mat$first, mat$second), alpha, conserv),
      updateC(
        mat$first,
        mat$second,
        f4sum,
        f2,
        f3,
        cross,
        alpha,
        conserv),
      ...
    )
  }


test_that("update works correctly", {
  compare_updates(1:100, FALSE, TRUE)

  A <- initializeR(cor(X, Y[, 192]), 0.05, FALSE)
  compare_updates(A, FALSE, TRUE)

  A <- 1:dx
  compare_updates(A, FALSE, TRUE)
})
