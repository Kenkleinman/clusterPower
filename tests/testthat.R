library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

#--------------------------------- SIMPLE CONTINUOUS OUTCOMES

context("Simple, normal outcome accuracy")

# compare to a reference value
test_that("continuous case matches reference", {
  expect_equal(round(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             nclusters = 5, nsubjects = 150,
                                             icc = .05, vart = 1, method = "weighted", 
                                             tol = .Machine$double.eps^0.25)), 7), as.numeric(0.4804988))
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
# NOTE: I'm not clear on the scale of the "difference" param vs. d
# and it was not clear from the code --- ASK KEN ABOUT THIS

test_that("simulation and analytic methods give similar power estimations", {
  sim.power <- cps.normal(nsim = 100, nsubjects = 150, nclusters = 5, difference = 4.804988,
                       ICC = 0.05, sigma = 100, alpha = 0.05, method = 'glmm', 
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


