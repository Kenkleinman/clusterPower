library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

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

#####################   FIXME
# needs to be parameterized properly
# analytic method power est is not inside conf.int. produced by simulation
# I still don't fully understand how "d" and "difference" are different measures of est treatment effect

test_that("simulation and analytic methods give similar power estimations", {
  sim.power <- cps.normal(nsim = 1000, nsubjects = 150, nclusters = 5, difference = 0.4804988,
                       ICC = 0.05, sigma = 100, alpha = 0.05, method = 'glmm', 
                       quiet = FALSE, all.sim.data = FALSE)
  reference <- as.numeric(round(crtpwr.2mean(alpha = 0.05, power = NA, d = 0.4804988, 
                            nclusters = 5, nsubjects = 150,
                            icc = .05, vart = 1, method = "weighted", 
                            tol = .Machine$double.eps^0.25), 1))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})

###################### END FIXME

context("Simple, binary outcome accuracy")

# does not match external value FIXME

test_that("binary case is accurate", {
  expect_equal(ceiling(as.numeric(crtpwr.2prop(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 150, cv = 0, 
                                             p1 = .1, p2 = .2, icc = .05, pooled = FALSE, 
                                             p1inc = TRUE, tol = .Machine$double.eps^0.25))), 
               ceiling(n4props(pe=.2, pc=.1, m=150, ICC=.05)$n))
})


