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
test_that("analytic normal irgtt case matches simulated method normal irgtt case", {
  expect_equal(signif(as.numeric(cps.irgtt.normal(nsim = 100, 
                                           nsubjects = c(400, 10),
                                           nclusters = 30, 
                                           difference = 0.5,
                                           sigma_sq = 4, 
                                           sigma_sq2 = 1, 
                                           sigma_b_sq2 = 0.5, 
                                           alpha = 0.05,
                                           quiet = FALSE, 
                                           all.sim.data = FALSE,
                                           seed = 123)$power$Power), 1), 
               as.numeric(signif(cpa.irgtt.normal(nclusters = 30, 
                                           nsubjects = 10, 
                                           d = 0.5, 
                                           varu = 0.5, 
                                           vare = 1, 
                                           varr = 4, 
                                           ncontrols = 400, 
                                           power = NA), 1)))})
