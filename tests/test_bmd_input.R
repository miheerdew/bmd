library(Rcpp)
library(testthat)
library(microbenchmark)

setwd("~/Dev/Bmd")
sourceCpp("bmd_input.cpp")
source("varcalcs.R")
source("sim_eQTL_network.R")

context("BmdInput")

sim <-sim_eQTL_network(make_param_list(b=10))
X <- scale(sim$X)
Y <- scale(sim$Y)
dy <- ncol(Y)
dx <- ncol(X)
n <- nrow(X)
inp <- new(BmdInput, X, Y)
alpha = 0.05

#test_that <- function(...) NA

test_that("BmdInput object meta data", {
  expect_equal(inp$n, nrow(X))
  expect_equal(inp$dx, ncol(X))
  expect_equal(inp$dy, ncol(Y))
  expect_equivalent(inp$X,X)
  expect_equivalent(inp$Y,Y)
 })

test_that("Computes cross_cors correctly", {
   expect_equal(inp$cross_cors(3:5, 5:10), cor(sim$X[,3:5],sim$Y[,5:10]))
   expect_equal(inp$cross_cors(5, 2), cor(sim$X[,5,drop=FALSE],sim$Y[,2,drop=FALSE]))
})

compare_variance <- function(A, test_x, benchmark=FALSE,...){
  if(test_x){
    fixdMat <- sim$Y[,A]
    thisMat <- sim$X
  } else {
    fixdMat <- sim$X[,A]
    thisMat <- sim$Y
  }
  if(benchmark){
    microbenchmark(C=inp$vars(A,test_x),
    R=varcalc1_multi(thisMat, fixdMat),...)
  } else {
    v1 <- inp$vars(A,test_x)
    v2 <- varcalc1_multi(thisMat, fixdMat)
    expect_equivalent(v1,v2)
    list(C=v1, R=v2, cormat=cor(thisMat, fixdMat))
  }
}

test_that("Variance is non-negative", {
  expect_true(all(inp$vars(1:3,TRUE)>=0))
  expect_true(all(inp$vars(1:20,FALSE)>=0))
})

test_that("Computes variance correctly", {
  compare_variance(1:3, FALSE)
  compare_variance(1:3, TRUE)
  compare_variance(1, FALSE)
  compare_variance(4, TRUE)
})

which_side <- function (A, test_x) {
  if (test_x) {
    return(list(fixdMat=sim$Y[,A],
            thisMat=sim$X))
  } else {
    return(list(fixdMat=sim$X[,A],
                thisMat=sim$Y))
  }
}
compare_pvalues <- function(A, test_x, benchmark=FALSE, ...){
 pval1 <- function() {
    mat <- which_side(A, test_x)
    vars <- varcalc1_multi(mat$thisMat, mat$fixdMat)
    zstats <- sqrt(n) * rowSums(cor(mat$thisMat, mat$fixdMat)) / sqrt(vars)
    pval <- as.vector(pnorm(zstats, lower.tail = FALSE))
 }

 if(benchmark){
   microbenchmark(R=pval1(), C=inp$pvals(A, test_x),...)
 } else {
  expect_equivalent(pval1(), as.vector(inp$pvals(A, test_x)))
 }
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

test_that("BH reject works correctly", {
  pvals <- runif(10000)
  p1 <- as.vector(bh_rejectC(pvals, 0.05, FALSE))
  p2 <- bh_rejectR(pvals, 0.05, FALSE)
  expect_equivalent(p1, p2)
  pvals <- c(runif(10000), 0, 1, 0)
  expect_equivalent(as.vector(bh_rejectC(pvals, 0.05, TRUE)),bh_rejectR(pvals, 0.05, TRUE))

  pvals <- as.vector(cor(X, Y[,23]))
  expect_equivalent(as.vector(bh_rejectC(pvals, 0.05, FALSE)),bh_rejectR(pvals, 0.05, FALSE))
  expect_equivalent(as.vector(bh_rejectC(pvals, 0.05, TRUE)),bh_rejectR(pvals, 0.05, TRUE))

  })

initialize3 <- function(u, test_x, conserv=TRUE){
  if(test_x){
      cor_to_u <- cor(X, Y[,u])
  } else {
      cor_to_u <- cor(Y, X[,u])
  }
  fischer_tranformed_cor <- atanh(cor_to_u)*sqrt(n-3)
  pvals <- pnorm(fischer_tranformed_cor, lower.tail = FALSE)
  return(bh_rejectR(pvals, alpha, conserv))
}

test_that("initialize works correctly", {
  res <- as.vector(inp$init(23,TRUE,TRUE))
  expect_equivalent(initialize3(23,TRUE,TRUE), res)

  res <- as.vector(inp$init(52,FALSE,TRUE))
  expect_equivalent(initialize3(52,FALSE,TRUE), res)

  res <- as.vector(inp$init(83, FALSE ,FALSE))
  expect_equivalent(initialize3(83, FALSE,FALSE), res)

  res <- as.vector(inp$init(122,TRUE,FALSE))
  expect_equivalent(initialize3(122,TRUE, FALSE), res)
})

pvals <- function(A, test_x) {
  if (test_x) {
    corsums <- as.vector(rowSums(inp$cross_cors(1:dx, A)))
  } else {
    corsums <- as.vector(colSums(inp$cross_cors(A, 1:dy)))
  }
  vars <- inp$vars(A, test_x)
  zstats <- sqrt(n) * corsums / sqrt(vars)
  return(pnorm(zstats, lower.tail = FALSE))
}

update5 <- function (B, test_x, conserv) {
  return(bh_rejectR(pvals(B, test_x), alpha, conserv))

}

test_that("update works correctly", {
  A <- inp$init(23,TRUE,TRUE)
  res <- as.vector(inp$update(A,FALSE,TRUE))
  expect_equivalent(update5(A,FALSE,TRUE), res)

  A <- inp$init(192,FALSE,TRUE)
  res <- as.vector(inp$update(A,TRUE,TRUE))
  expect_equivalent(update5(A,TRUE,TRUE), res)

  A <- 1:dx
  res <- as.vector(inp$update(A,FALSE,TRUE))
  expect_equivalent(update5(A,FALSE,TRUE), res)
})



pvals <- runif(10000)
res <- microbenchmark(  pC = as.vector(bh_rejectC(pvals, 0.05, FALSE)),
                        pR = bh_rejectR(pvals, 0.05, FALSE), times=20)
print(res)
