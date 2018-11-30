library(testthat)
library(clusterPower)

test_check("clusterPower")

context("Simple, normal outcome accuracy")

# compare to a reference value
test_that("continuous case is accurate", {
  expect_equal(round(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             nclusters = 5, nsubjects = 150,
                                             icc = .05, vart = 1, method = "weighted", 
                                             tol = .Machine$double.eps^0.25)), 7), as.numeric(0.4804988))
})

#####################   FIXME
# needs to be parameterized properly
# analytic method power est is not inside conf.int. produced by simulation
# I still don't fully understand how "d" and "difference" are different measures of est treatment effect

test_that("simulation and analytic methods give similar power estimations"), {
  holder <- cps.normal(nsim = 100, nsubjects = 150, nclusters = 5, difference = .1857,
                       ICC = 0.05, sigma = 100, alpha = 0.05, method = 'glmm', 
                       quiet = FALSE, all.sim.data = FALSE)
  expect_equal(as.numeric(round(crtpwr.2mean(alpha = 0.05, power = NA, d = .285, 
                            nclusters = 5, nsubjects = 150,
                            icc = .05, vart = 1, method = "weighted", 
                            tol = .Machine$double.eps^0.25), 1)), as.numeric(holder$power$Power))
}
###################### END FIXME

context("Simple, binary outcome accuracy")

# compare to a reference value -- should be replaced with an external value
test_that("binary case is accurate", {
  expect_equal(round(as.numeric(crtpwr.2prop(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 150, cv = 0, 
                                             p1 = .1, p2 = .2, icc = .05, pooled = FALSE, 
                                             p1inc = TRUE, tol = .Machine$double.eps^0.25)), 5), as.numeric(11.05385))
})

