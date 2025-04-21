library("Matrix")
library("testthat")
library("L0TFinv")


test_that("The sample size is less than the matrix order , making the problem underdetermined", {
  expect_error(L0TFinv::DiffMat(n = 2, q = 2))

  expect_error(L0TFinv::L0TFinv.opt(y=c(1,2,3,5), kmax=10, q=0))
  expect_error(L0TFinv::L0TFinv.fix(y=c(-10,10), k=3, q=1))

  expect_error(L0TFinv::SimuBlocksInv(n = 2, sigma = 0.1, tau = c(0.1,0.9), h = c(1,2,3)))
  expect_error(L0TFinv::SimuWaveInv(n = 2, sigma = 0.1, tau = c(0.1,0.9), h = c(1,2,3), a0 = -10))
})


test_that("Active set exceeding the upper/lower bounds", {
  A0 <- c(0,3,10)
  B0 <- c(3,10,15,20)
  expect_error(L0TFinv::solMat(n = 15, q = 0, A = A0))
  expect_error(L0TFinv::solMat(n = 15, q = 1, A = B0))

  expect_error(L0TFinv::SimuBlocksInv(n = 100, sigma = 0.1, tau = c(-0.1,0.9), h = c(1,2,3)))
  expect_error(L0TFinv::SimuWaveInv(n = 100, sigma = 0.1, tau = c(-0.1,0.9), h = c(1,2,3), a0 = -10))
})


test_that("Matrix order is neither 0 nor 1", {
  A1 = c(1,2,5,8)
  A2 = c(1,3,8,10,15)
  expect_error(L0TFinv::solMat(n = 10, q = 3, A = A1))
  expect_error(L0TFinv::solMat(n = 15, q = 5, A = A2))
})


## Test sMatrix.R
test_that("Generate specialized matrix", {
  Mat1 <- L0TFinv::DiffMat(n = 10, q = 0)
  Mat2 <- L0TFinv::DiffMat(n = 15, q = 1)
  Mat3 <- L0TFinv::DiffMat(n = 15, q = 2)
  expect_equal(dim(Mat1), c(9, 10))
  expect_equal(dim(Mat2), c(13, 15))
  expect_equal(dim(Mat3), c(12, 15))

  Xmat1 <- L0TFinv::XMat(n = 15, q = 0)
  Xmat2 <- L0TFinv::XMat(n = 20, q = 1)
  Xmat3 <- L0TFinv::XMat(n = 25, q = 2)
  expect_equal(dim(Xmat1), c(15, 15))
  expect_equal(dim(Xmat2), c(20, 20))
  expect_equal(dim(Xmat3), c(25, 25))
})


## Test solveMatrix.R
test_that("Compute the inversion of the specialized matrix", {
  A1 = c(1,2,5,8)
  S1 <- L0TFinv::solMat(n = 10, q = 0, A = A1)
  A2 = c(1,3,8,10,15)
  S2 <- L0TFinv::solMat(n = 15, q = 1, A = A2)
  expect_equal(dim(S1), c(4, 4))
  expect_equal(dim(S2), c(5, 5))

  S3 <- L0TFinv::solMat(n = 10, q = 0, A = c(3))
  S4 <- L0TFinv::solMat(n = 15, q = 0, A = c(4,8))
})


## Test invL0TF.R
test_that("Model fitting for piecewise constant trends", {
  tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
  h = c(-1, 5, 3, 0, -1, 2)
  n = 500
  BlocksData <- L0TFinv::SimuBlocksInv(n = n, sigma = 0.2, seed = 50, tau = tau ,h = h)
  resoptbic <- L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=20, q=0, first=0.01, last=1, penalty="bic")
  resoptsic <- L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=20, q=0, first=0.01, last=1, penalty="sic")
  resfix <- L0TFinv::L0TFinv.fix(y=BlocksData$y, k=6, q=0, first=0, last=1)
  expect_no_error(print(resoptbic))
  expect_no_error(print(resoptsic))
  expect_no_error(print(resfix))
  expect_equal(length(resfix$Ak),6)

  expect_error(L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=0.01, last=1, penalty="log(n)"))
  expect_error(L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=0.01, last=2))
  expect_error(L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=-0.01, last=1))
  expect_error(L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=10, q=2, first=0, last=1))

  expect_error(L0TFinv::L0TFinv.fix(y=BlocksData$y, k=6, q=0, first=0, last=2))
  expect_error(L0TFinv::L0TFinv.fix(y=BlocksData$y, k=6, q=0, first=-1, last=1))
  expect_error(L0TFinv::L0TFinv.fix(y=BlocksData$y, k=6, q=2, first=0, last=1))

  expect_error(plot(resfix,type="all"))
  expect_no_error(plot(resfix,type="mse"))
  expect_no_error(plot(resfix,type="bic"))
  expect_no_error(plot(resfix,type="sic"))
  expect_no_error(plot(resfix,type="yhat",k=5))
  expect_no_error(plot(resfix,type="yhat"))
  expect_error(plot(resfix,type="yhat",k=7))
})


test_that("Model fitting for piecewise linear trends", {
  tau1 = c(0.4, 0.6, 0.7)
  h1 = c(-3, 5, -4, 6)
  a0 = -10
  WaveData <- L0TFinv::SimuWaveInv(n = 500, sigma = 0.1, seed = 50, tau = tau1, h = h1, a0 = a0)
  res1 <- L0TFinv::L0TFinv.opt(y=WaveData$y, kmax=10, q=1, first=0, last=0.99, penalty="sic")
  res1k <- L0TFinv::L0TFinv.fix(y=WaveData$y, k=8, q=1, first=0, last=0.99)
  expect_no_error(coef(res1))
  expect_no_error(coef(res1,k=7))
  expect_error(coef(res1,k=11))
  expect_no_error(coef(res1k))
  expect_no_error(coef(res1k,k=6))
  expect_error(coef(res1k,k=9))
  expect_equal(length(res1k$Ak),8)

  expect_error(plot(res1,type="all"))
  expect_no_error(plot(res1,type="mse"))
  expect_no_error(plot(res1,type="bic"))
  expect_no_error(plot(res1,type="sic"))
  expect_no_error(plot(res1,type="yhat",k=5))
  expect_no_error(plot(res1,type="yhat"))
  expect_error(plot(resfix,type="yhat",k=12))
})


## Test TFmetrics.R
test_that("Generate the metrics for the trend filtering model", {
  tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
  h = c(-1, 5, 3, 0, -1, 2)
  n = 500
  BlocksData <- L0TFinv::SimuBlocksInv(n = n, sigma = 0.2, seed = 50, tau = tau ,h = h)
  res <- L0TFinv::L0TFinv.opt(y=BlocksData$y, kmax=10, q=0, first=0.01, last=1, penalty="bic")
  metrics <- L0TFinv::TFmetrics(BlocksData$y0,BlocksData$tau,res$yopt,res$Aopt/n)
  metrics0 <- L0TFinv::TFmetrics(y0 = BlocksData$y0,tau = BlocksData$tau,yhat = res$yopt)

  expect_equal(length(metrics), 4)
  expect_equal(length(metrics0), 2)
})


## Test dataSimu.R
test_that("Generate the metrics for the trend filtering model", {
  tau = c(0.1, 0.3, 0.4, 0.7, 0.85)
  h = c(-1, 5, 3, 0, -1, 2)
  tau1 = c(0.4, 0.6, 0.7)
  h1 = c(-3, 5, -4, 6)
  a0 = -10
  n = 500
  BlocksData <- L0TFinv::SimuBlocksInv(n = n, sigma = 0.2, tau = tau ,h = h)
  WaveData <- L0TFinv::SimuWaveInv(n = n, sigma = 0.1, tau = tau1, h = h1, a0 = a0)
  expect_equal(length(BlocksData$tau), 5)
  expect_equal(length(WaveData$tau), 3)

  expect_error(L0TFinv::SimuBlocksInv(n = 100, sigma = 0.1, tau = c(0.1,0.9), h = c(1,4)))
  expect_error(L0TFinv::SimuWaveInv(n = 100, sigma = 0.1, tau = c(-0.1,0.9), h = c(2,3), a0 = -10))
})


