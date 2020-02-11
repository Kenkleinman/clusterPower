library(testthat)
library(clusterPower)
library(CRTSize)

library("devtools")
install_github("Kenkleinman/clusterPower@lexi")

library("clusterPower")

test_check("clusterPower")

#validateVariance 
#optimizerSearch


############################################
######  TESTING FOR SMALL FUNCTIONS  #######
############################################

context("confint.calc")

test_that("confint.calc matches a reference", {
  some_nums <- confint.calc(nsim = 1000, alpha = 0.05, p.val = 0.01, 
                            names.power = c("bestCI"))
  expect_equal(3, length(some_nums))
  expect_equal(0.001, some_nums$Power)
  expect_equal(2.531749e-05, some_nums$Lower.95.CI)
  expect_equal(0.005558924, some_nums$Upper.95.CI)
})

context("createMissingVarianceParam")

test_that("createMissingVarianceParam matches a reference", {
  ICC <- createMissingVarianceParam(sigma_sq = c(1, 1, 0.9), sigma_b_sq = c(0.1, 0.15, 0.1))
  expect_equal(3, length(ICC))
  expect_equal(0.09090909, ICC[1])
  expect_equal(0.13043478, ICC[2])
  expect_equal(0.10000000, ICC[3])
  sig_b <- createMissingVarianceParam(sigma_sq = c(1, 1, 0.9), ICC = ICC)
  expect_equal(3, length(sig_b))
  expect_equal(0.10, sig_b[1])
  expect_equal(0.15, sig_b[2])
  expect_equal(0.10, sig_b[3])
})

context("is.wholenumber")

test_that("is.wholenumber matches a reference", {
  expect_equal(FALSE, is.wholenumber(3.5))
  expect_equal(TRUE, is.wholenumber(5001))
})


context("prop_H0_rejection")

test_that("prop_H0_rejection matches a reference", {
  prop <- prop_H0_rejection(alpha = 0.05, nsim = 1000, 
                    LRT.holder.abbrev = 804)
  expect_equal(3, length(prop))
  expect_equal(0.804, as.numeric(prop[1]))
  expect_equal(0.828, as.numeric(prop[3]))
  expect_equal("Power", names(prop[1]))
})


# this is still messed up
context("type1ErrTest")

test_that("type1ErrTest matches a reference", {
 type1ErrTest(sigma_sq_ = c(0.1, 4),
                       sigma_b_sq_ = c(0.1, 0.15),
                       nsubjects_ = list(rep(4, 10),
                                         rep(4, 10)))
                       
                       
  expect_equal(3, length(prop))
  expect_equal(0.804, as.numeric(prop[1]))
  expect_equal(0.828, as.numeric(prop[3]))
  expect_equal("Power", names(prop[1]))
})

if ((nobstemp < 5 & ICCtemp < 0.1 & clusterstemp < 20) | 
    (nobstemp < 10  & ICCtemp < 0.05 & clusterstemp < 20) |
    (nobstemp < 20  & ICCtemp < 0.02 & clusterstemp < 20) |
    (nobstemp < 50  & ICCtemp < 0.01 & clusterstemp < 10))
  
  createMissingVarianceParam(sigma_sq = c(0.1, 4), sigma_b_sq = c(0.001, 0.15))
