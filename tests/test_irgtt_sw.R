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
  expect_equal(as.numeric(cps.irgtt.normal(nsim = 100, 
                                           nsubjects = 10,
                                           nclusters = 10, 
                                           difference = .5,
                                           sigma_sq = 1, 
                                           sigma_sq2 = 0.9, 
                                           sigma_b_sq2 = 0.1, 
                                           alpha = 0.05,
                                           quiet = FALSE, 
                                           all.sim.data = FALSE)), 
               as.numeric(cpa.irgtt.normal(nclusters = 10, 
                                           nsubjects = 10, 
                                           d = 0.5, 
                                           varu = 0.1, 
                                           vare = 0.9, 
                                           varr = 1, 
                                           ncontrols = 77.81, 
                                           power = NA)))})
