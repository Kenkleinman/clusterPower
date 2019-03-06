library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

#######################################
#### MULTI-ARM CRT DESIGNS TESTING #######
#######################################

#--------------------------------- MULTI-ARM CONTINUOUS OUTCOMES

context("Multi-arm, normal outcome accuracy")

# compare to a reference value from NIH calculator
test_that("binary multi-arm case matches 2-arm binary case (simulated method)", {
  expect_equal(as.numeric(cps.ma.binary(nsim = 100, nsubjects = 50, narms=2,
                                                nclusters=4,
                                                probs = c(0.4, 0.1),
                                                sigma_b_sq = 1, alpha = 0.05,
                                                quiet = FALSE, method = 'glmm', 
                                                all.sim.data = FALSE, 
                                                multi.p.method="none",
                                                poor.fit.override = TRUE,
                                                seed = 123, 
                                                cores="all")[[1]][1]), 
                       as.numeric(cps.binary(nsim = 100, nsubjects = 50, 
                                             nclusters = 4, p1 = 0.4,p2 = 0.1, 
                                             sigma_b = 1, alpha = 0.05, method = 'glmm', 
                                             all.sim.data = FALSE, seed = 123)[[3]][1]))})

#FIXME here's where I left off.
# compare to a comparable function in CRTSize -uses ICC
test_that("continuous case matches CRTSize", {
  expect_equal(round(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 150,
                                             icc = .05, vart = 1, method = "weighted", 
                                             tol = .Machine$double.eps^0.25, d=0.4804988)), 0), 
               round(n4means(delta=0.4804988, sigma=1, m=150, ICC=0.05)$n),0)
})

# compare simulation and analytic methods for continuous outcomes -uses ICC, asked John in include the definition in the documentation
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

