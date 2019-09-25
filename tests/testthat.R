library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

#######################################
#### SIMPLE CRT DESIGNS TESTING #######
#######################################

#--------------------------------- SIMPLE CONTINUOUS OUTCOMES

context("Simple, normal outcome accuracy")

# compare to a reference value from NIH calculator
test_that("continuous case matches reference", {
  expect_equal(ceiling(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             d=0.4804988, nsubjects = 150,
                                             icc = .05, vart = 1, method = "weighted", 
                                             tol = .Machine$double.eps^0.25))), as.numeric(5))
})

# compare to a comparable function in CRTSize
test_that("continuous case matches CRTSize", {
  expect_equal(round(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 150,
                                             icc = .05, vart = 1, method = "weighted", 
                                             tol = .Machine$double.eps^0.25, d=0.4804988)), 0), 
               round(n4means(delta=0.4804988, sigma=1, m=150, ICC=0.05)$n),0)
})

# compare simulation and analytic methods for continuous outcomes
test_that("simulation and analytic methods give similar power estimations", {
  sim.power <- cps.normal(nsim = 100, nsubjects = 150, nclusters = 5, difference = .4804988,
                       ICC = 0.05, sigma = 1, alpha = 0.05, method = 'glmm', 
                       quiet = FALSE, all.sim.data = FALSE)
  print(sim.power$power)
  reference <- as.numeric(round(crtpwr.2mean(alpha = 0.05, power = NA, d = 0.4804988, 
                            nclusters = 5, nsubjects = 150,
                            icc = .05, vart = 1, method = "weighted", 
                            tol = .Machine$double.eps^0.25), 1))
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})



#--------------------------------- SIMPLE BINARY OUTCOMES

context("Simple, binary outcome accuracy")

# NOTE: does not match NIH reference value. Seems to consistently underestimate the nclusters needed 
# by 1. NIH estimates nclusters using a t-distribution. Do we need them to match more closely?

# Ken says add user option to use t or normal distribution.

test_that("binary case matches reference value from NIH calculator", {
  expect_equal(ceiling(as.numeric(crtpwr.2prop(alpha = 0.05, power = 0.8, 
                                               nclusters = NA, nsubjects = 150, cv = 0, 
                                               p1 = .1, p2 = .2, icc = .05))), 13-1)
})

# p1=treatment, p2=control
# NOTE: does not match n4props value. Seems to consistently underestimate the nclusters needed 
# by 1. n4props estimates nclusters using a t-distribution. Do we need them to match more closely?

test_that("binary case matches CRTSize::n4props", {
  expect_equal(ceiling(as.numeric(crtpwr.2prop(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 150, cv = 0, 
                                             p1 = .1, p2 = .2, icc = .05))), 
               ceiling(n4props(pe=.2, pc=.1, m=150, ICC=.05)$n-1))
})

# compare simulation and analytic methods for binary outcomes
# NOTE: sigma_b has no default and may not be zero (although the help page
# says "if sigma_b2 is not specified, between cluster variances are
# assumed to be equal for both groups"). I'm probably not 
# understanding something. --- ASK KEN ABOUT THIS

# Ken says: Should be a test for sigma b2=NA and if Sigb2=NA val set to sigma b.

test_that("simulation and analytic methods give similar power estimations for binary outcomes", {
  sim.power <- cps.binary(nsim = 100, nsubjects = 150, nclusters = 12, p1 = 0.1,
                          p2 = 0.2, sigma_b = .4, alpha = 0.05, method = 'glmm', 
                          all.sim.data = FALSE)
  print(sim.power$power)
  reference <- round(as.numeric(crtpwr.2prop(alpha = 0.05, power = NA, 
                                               nclusters = 12, nsubjects = 150, cv = 0, 
                                               p1 = .1, p2 = .2, icc = .05)), 1)
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})


#--------------------------------- SIMPLE POISSON OUTCOMES

context("Simple, poisson (incidence rate) outcome accuracy")

# incidence rate comparison
test_that("incidence rate outcomes matches CRTSize::n4incidence", {
  expect_equal(ceiling(as.numeric(crtpwr.2rate(alpha = 0.05, power = 0.8, 
                                               nclusters = NA, py = 1000, 
                                               r1 = 0.01, r2 = 0.005, cvb = 0.25))), 
               ceiling(n4incidence(le=0.01, lc=0.005, m=1000, t=1, CV=0.25)$n))
})

# compare simulation and analytic methods for poisson outcomes
# test .7-.9 in increments of .25
# NOTE: c1 & c2 (treatment) are outcome counts while r1 (treatment) & r2 are outcome rates
test_that("simulation and analytic methods give similar power estimations for poisson outcomes", {
  sim.power <- cps.count(nsim = 100, nsubjects = 120, nclusters = 4, c1 = 20,
                         c2 = 10, sigma_b = 0.1, family = 'poisson', analysis = 'poisson',
                         method = 'glmm', alpha = 0.05, quiet = FALSE, all.sim.data = TRUE)
  print(sim.power$power)
  reference <- as.numeric(crtpwr.2rate(alpha = 0.05, power = NA, 
                                               nclusters = 4, py = 120, 
                                               r1 = 0.1, r2 = 0.2, cvb = 0.1))
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})


#######################################
###### DID CRT DESIGNS TESTING ########
#######################################

#--------------------------------- DID CONTINUOUS OUTCOMES

context("DID normal outcome accuracy")

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cross-sectional)", {
  expect_equal(ceiling(as.numeric(crtpwr.2meanD(alpha = 0.05, power = 0.8, 
                                               d = 0.48, nsubjects = 150,
                                               rho_c = 0, rho_s = 0,
                                               icc = .05, vart = 1))), 
               as.numeric(9))
})

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(crtpwr.2meanD(alpha = 0.05, power = 0.8, 
                                                d = 0.48, nsubjects = 150,
                                                rho_c = 0.3, rho_s = 0.7,
                                                icc = .05, vart = 1))), 
               as.numeric(7))
})


test_that("DID normal simulation and analytic methods give similar power estimations", {
  sim.power <- cps.did.normal(nsim = 100, nsubjects = 150, nclusters = 6, 
                              difference = .48, sigma = 1, alpha = 0.05, sigma_b0 = .1,
                              method = 'glmm', quiet = FALSE, all.sim.data = FALSE)
  print(sim.power$power)
  reference <- as.numeric(round(crtpwr.2meanD(alpha = 0.05, power = NA, d = 0.48, 
                                             nclusters = 6, nsubjects = 150,
                                             icc = .05, vart = 1, rho_c = 0.3, 
                                             rho_s = 0.7), 1))
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})


#--------------------------------- DID BINARY OUTCOMES

context("DID binary outcome accuracy")

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(crtpwr.2propD(nsubjects=50, p=.5, 
                                                d=.1, icc=.05, rho_c=.3, 
                                                rho_s=.7))), 
               as.numeric(33))
})

# NOTE: cps.did.binary seems to give unusually low estimates of power. I don't know why yet.
# NOTE: the simulation method is very slow
test_that("DID binary simulation and analytic methods give similar power estimations", {
  sim.power <- cps.did.binary(nsim = 10, nsubjects = 150, nclusters = 30, 
                              p1 = 0.4, p2 = 0.5, sigma_b0 = 1, alpha = 0.05,
                              method = 'glmm', all.sim.data = FALSE)
  print(sim.power$power)
  reference <- as.numeric(round(crtpwr.2propD(nsubjects=150, p=.4,
                                              power=NA,
                                              nclusters = 30, 
                                              d=.1, icc=.05, rho_c=0, 
                                              rho_s=0), 1))
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})


# still need to test:
# DID poisson outcomes
# IRGTT normal/binary/count
# SW normal/binary/count
# Multi-arm normal/binary/count




###########################################
######## VALIDATION FXNS TESTS ############
###########################################

testthat::context("createMissingVarianceParam: calculate ICC, sigma, 
                  or sigma_b based on user providing 2 of those values")

# compare to a reference value 
testthat::test_that("createMissingVarianceParam returns ICC", {
  testthat::expect_equal(createMissingVarianceParam(sigma = c(1, 1, 0.9), 
                                                       sigma_b = c(0.1, 0.15, 0.1)), 
                            c(0.09090909, 0.13043478, 0.10000000))})

testthat::test_that("createMissingVarianceParam fails when fewer than 2 params provided", {
  testthat::show_failure(createMissingVarianceParam(sigma = c(1, 1, 0.9)))})

  