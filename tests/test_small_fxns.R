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


context("type1ErrTest")

test_that("type1ErrTest matches a reference", {
 warn <- type1ErrTest(sigma_sq_ = c(0.1, 4),
                       sigma_b_sq_ = c(0.1, 0.15),
                       nsubjects_ = list(rep(4, 10),
                                         rep(4, 10)))
 expect_equal(310, nchar(warn))
 nowarn <- type1ErrTest(sigma_sq_ = c(1, 1),
                      sigma_b_sq_ = c(0.5, 0.1),
                      nsubjects_ = list(rep(20, 20),
                                        rep(20, 20)))
 expect_equal(TRUE, is.null(nowarn))
 })


context("validateVariance")

test_that("validateVariance matches a reference", {
  validateVariance(dist = NULL, difference = NULL, alpha = alpha, 
   ICC = ICC, sigma_sq = sigma_sq, 
   sigma_b_sq = sigma_b_sq, ICC2 = NA, sigma_sq2 = NA, 
   sigma_b_sq2 = NA, method = method, quiet = quiet, 
   all.sim.data = all.sim.data, 
   poor.fit.override = poor.fit.override, 
   cores = NA,
   probs = NA)
  expect_equal(TRUE, is.null(nowarn))
})
