library(testthat)
library(clusterPower)
library(CRTSize)

library("devtools")
install_github("Kenkleinman/clusterPower@lexi")

library("clusterPower")

test_check("clusterPower")



#######################################
######  IRGTT METHODS TESTING   #######
#######################################

#--------------------------------- IRGTT NORMAL OUTCOME

context("IRGTT outcome accuracy")

# compare analytic to simulation methods
test_that("analytic normal irgtt case matches simulated method irgtt case", {
  expect_equal(as.numeric(cps.ma.binary(nsim = 100, nsubjects = 50, narms = 2,
                                        nclusters = 4,
                                        probs = c(0.4, 0.1),
                                        sigma_b_sq = 1, alpha = 0.05,
                                        quiet = FALSE, method = 'glmm', 
                                        all.sim.data = FALSE, 
                                        multi.p.method = "none",
                                        poor.fit.override = TRUE,
                                        seed = 123, 
                                        cores = "all")[[1]][1]), 
               as.numeric(crtpwr.irgtt.mean(nclusters=10, 
                                            nsubjects = 10, 
                                            d = 0.5, 
                                            varu = 0.1, 
                                            vare = 0.9, 
                                            varr = 1)))})
