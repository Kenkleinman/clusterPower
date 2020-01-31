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

# compare analytic to simulation methods (normal)
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


# compare simulation method (normal) to a constant
test_that("simulation normal irgtt case matches a constant", {
  expect_equal(signif(as.numeric(cps.irgtt.normal(nsim = 100, 
                                                  nsubjects = c(550, 12),
                                                  nclusters = 125, 
                                                  difference = 0.25,
                                                  sigma_sq = 3, 
                                                  sigma_sq2 = 2.1, 
                                                  sigma_b_sq2 = 0.25, 
                                                  alpha = 0.05,
                                                  quiet = FALSE, 
                                                  all.sim.data = FALSE,
                                                  seed = 123)$power$Power), 2),
               0.78)})

# compare analytic to a constant (normal)
test_that("analytic normal irgtt case matches a constant", {
  expect_equal(as.numeric(signif(cpa.irgtt.normal(nclusters = 125, 
                                                  nsubjects = 12, 
                                                  d = 0.25, 
                                                  varu = 0.25, 
                                                  vare = 2.1, 
                                                  varr = 3, 
                                                  ncontrols = 550, 
                                                  power = NA), 2)),
               0.76)})

#--------------------------------- IRGTT BINARY OUTCOME

# compare analytic to simulation methods (binary)
test_that("analytic binary irgtt case matches simulated method binary irgtt case", {
  expect_equal(signif(as.numeric(cps.irgtt.binary(nsim = 100, 
                                                  nsubjects = 30, 
                                                  nclusters = 10, 
                                                  p1 = 0.2,
                                                  p2 = 0.5, 
                                                  sigma_b_sq2 = 1, 
                                                  alpha = 0.05, 
                                                  all.sim.data = FALSE, 
                                                  seed = 123)$power$Power), 1), 
               as.numeric(signif(cpa.irgtt.binary(nclusters = 10, 
                                                  nsubjects = 30,
                                                  ncontrols = 30, 
                                                  icc = 0.25, 
                                                  p.e = 0.5, 
                                                  p.c = 0.2, 
                                                  power = NA), 1)))})

# compare simulation methods (binary) to a constant
test_that("simulated method binary irgtt case matches a constant", {
  expect_equal(signif(as.numeric(cps.irgtt.binary(nsim = 100, 
                                                  nsubjects = c(300, 12), 
                                                  nclusters = 120, 
                                                  p1 = 0.3,
                                                  p2 = 0.2, 
                                                  sigma_b_sq2 = 1.75, 
                                                  alpha = 0.05, 
                                                  all.sim.data = FALSE, 
                                                  seed = 123)$power$Power), 2), 
               0.81)})


# compare analytic methods (binary) to a constant
test_that("analytic binary irgtt case matches a constant", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = 0.4, 
                                                  p.e = 0.2, 
                                                  p.c = 0.3, 
                                                  power = NA), 2)),
               0.79)})


###############################################
######  STEPPED-WEDGE METHODS TESTING   #######
###############################################

#--------------------------------- SW NORMAL OUTCOME

context("SW outcome accuracy")

# compare analytic to simulation methods (normal)
test_that("analytic normal SW case matches simulated method normal SW case", {
  expect_equal(signif(as.numeric(cps.sw.normal(nsim = 100, nsubjects = 50, nclusters = 30, 
                                               difference = 1.75, steps = 5, sigma_sq = 100, sigma_b_sq = 30, 
                                               alpha = 0.05, method = 'glmm', quiet = FALSE, 
                                               all.sim.data = FALSE)$power$Power), 1), 
               as.numeric(signif(
                 
                 cpa.sw.normal(nsubjects = 50, 
                               nclusters = 300, 
                               ntimes = 5, 
                               d = 1.75, 
                               icc = 0.5, 
                               rho_c = 100, 
                               rho_s = 30, 
                               vart = 130, 
                               power = NA)
                 
               ), 1)))})


