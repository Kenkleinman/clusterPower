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
    sim_model <- cps.irgtt.normal(nsim = 100, 
                                           nsubjects = c(400, 10),
                                           nclusters = 30, 
                                           difference = 0.5,
                                           sigma_sq = 4, 
                                           sigma_sq2 = 1, 
                                           sigma_b_sq2 = 0.5, 
                                           alpha = 0.05,
                                           quiet = FALSE, 
                                           all.sim.data = FALSE,
                                           seed = 123) 
               analytic_model <- cpa.irgtt.normal(nclusters = 30, 
                                           nsubjects = 10, 
                                           d = 0.5, 
                                           varu = 0.5, 
                                           vare = 1, 
                                           varr = 4, 
                                           ncontrols = 400, 
                                           power = NA)
                                 expect_equal(TRUE, data.table::between(analytic_model,
                                                                        sim_model$power$Lower.95.CI,
                                                                        sim_model$power$Upper.95.CI))
                                 })


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
test_that("analytic binary irgtt power is within simulated method binary irgtt CI", {
    sim_model <- cps.irgtt.binary(nsim = 100, 
                                                  nsubjects = 30, 
                                                  nclusters = 10, 
                                                  p1 = 0.2,
                                                  p2 = 0.5, 
                                                  sigma_b_sq2 = 1, 
                                                  alpha = 0.05, 
                                                  all.sim.data = FALSE, 
                                                  seed = 123)
               analytic_model <- cpa.irgtt.binary(nclusters = 10, 
                                                  nsubjects = 30,
                                                  ncontrols = 30, 
                                                  icc = 0.25, 
                                                  p.e = 0.5, 
                                                  p.c = 0.2, 
                                                  power = NA)
               expect_equal(TRUE, data.table::between(analytic_model, 
                                                      sim_model$power$Lower.95.CI,
                                                      sim_model$power$Upper.95.CI))
               })

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
test_that("analytic binary irgtt case matches a constant: power", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = 0.4, 
                                                  p.e = 0.2, 
                                                  p.c = 0.3, 
                                                  power = NA), 2)),
               0.79)})

test_that("analytic binary irgtt case matches a constant: icc calc", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = NA, 
                                                  p.e = 0.2, 
                                                  p.c = 0.3, 
                                                  power = 0.79), 2)),
               0.41)})

test_that("analytic binary irgtt case matches a constant: nsubjects", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = NA,
                                                  ncontrols = 300, 
                                                  icc = 0.4054797, 
                                                  p.e = 0.2, 
                                                  p.c = 0.3, 
                                                  power = 0.79), 2)),
               12)})

test_that("analytic binary irgtt case matches a constant: nclusters", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = NA, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = 0.4054797, 
                                                  p.e = 0.2, 
                                                  p.c = 0.3, 
                                                  power = 0.79), 2)),
               120)})

test_that("analytic binary irgtt case matches a constant: p.c", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = 0.4054797, 
                                                  p.e = 0.2, 
                                                  p.c = NA, 
                                                  power = 0.79), 2)),
               0.3)})

test_that("analytic binary irgtt case matches a constant: p.e", {
  expect_equal(as.numeric(signif(cpa.irgtt.binary(nclusters = 120, 
                                                  nsubjects = 12,
                                                  ncontrols = 300, 
                                                  icc = 0.4054797, 
                                                  p.e = NA, 
                                                  p.c = 0.3, 
                                                  power = 0.79), 2)),
               0.2)})

#--------------------------------- IRGTT COUNT OUTCOME

# compare simulation methods (count) to a constant
test_that("analytic irgtt power is within simulated method count irgtt CI", {

  irgtt.count.sim <- cps.irgtt.count(nsim = 100, nsubjects = c(500, 10), nclusters = 100, 
                                c1 = 85, c2 = 450, sigma_b_sq2 = 25, 
                                family = 'poisson', analysis = 'poisson',
                                alpha = 0.05, quiet = FALSE, all.sim.data = FALSE, 
                                opt = "L-BFGS-B")
  
    sim_model <- cps.irgtt.count(nsim = 100, 
                                                  nsubjects = 30, 
                                                  nclusters = 10, 
                                                  p1 = 0.2,
                                                  p2 = 0.5, 
                                                  sigma_b_sq2 = 1, 
                                                  alpha = 0.05, 
                                                  all.sim.data = FALSE, 
                                                  seed = 123) 
  analytic_model <- cpa.irgtt.binary(nclusters = 10, 
                                                  nsubjects = 30,
                                                  ncontrols = 30, 
                                                  icc = 0.25, 
                                                  p.e = 0.5, 
                                                  p.c = 0.2, 
                                                  power = NA)
expect_equal(TRUE, data.table::between(analytic_model, 
                                       sim_model$power$Lower.95.CI,
                                       sim_model$power$Upper.95.CI))
})


###############################################
######  STEPPED-WEDGE METHODS TESTING   #######
###############################################

#--------------------------------- SW NORMAL OUTCOME

context("SW outcome accuracy")

# compare analytic to simulation methods (normal)
test_that("analytic normal SW outcome is within estimated CI of simulated SW method", {
    sim_model <- cps.sw.normal(nsim = 100, 
                                               nsubjects = 50, 
                                               nclusters = 30, 
                                               difference = 1.75, 
                                               steps = 5, 
                                               sigma_sq = 100, 
                                               sigma_b_sq = 30, 
                                               alpha = 0.05, 
                                               method = 'glmm', 
                                               quiet = FALSE, 
                                               all.sim.data = FALSE)
               analytic_model <- cpa.sw.normal(nsubjects = 50, 
                                               nclusters = 300, 
                                               ntimes = 5, 
                                               d = 1.75, 
                                               icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                               rho_c = 0.5, 
                                               rho_s = 0.25, 
                                               vart = 2322.988, 
                                               power = NA)
               expect_equal(TRUE, 
                            data.table::between(analytic_model, 
                                                sim_model$power$Lower.95.CI, 
                                                sim_model$power$Upper.95.CI))
  })

# compare SW simulation method to a constant (normal)
test_that("simulated method normal SW case matches a constant", {
  expect_equal(signif(as.numeric(cps.sw.normal(nsim = 100, 
                                               nsubjects = 50, 
                                               nclusters = 30, 
                                               difference = 1.75, 
                                               steps = 5, 
                                               sigma_sq = 100, 
                                               sigma_b_sq = 30, 
                                               alpha = 0.05, 
                                               method = 'glmm', 
                                               quiet = FALSE, 
                                               all.sim.data = FALSE)$power$Power), 1), 
               0.8)})

# compare analytic SW methods (normal) to a constant
test_that("analytic normal SW case matches a constant", {
  expect_equal(as.numeric(signif(cpa.sw.normal(nsubjects = 50, 
                                               nclusters = 300, 
                                               ntimes = 5, 
                                               d = 1.75, 
                                               icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                               rho_c = 0.5, 
                                               rho_s = 0.25, 
                                               vart = 2322.988, 
                                               power = NA), 1)), 
               0.8)})

test_that("analytic normal SW case matches a constant: vart", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50, 
                                               nclusters = 300, 
                                               ntimes = 5, 
                                               d = 1.75, 
                                               icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                               rho_c = 0.5, 
                                               rho_s = 0.25, 
                                               vart = NA, 
                                               power = 0.8100001)), 
               2322.988)})

test_that("analytic normal SW case matches a constant: nsubjects", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = NA, 
                                        nclusters = 300, 
                                        ntimes = 5, 
                                        d = 1.75, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = 0.5, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               50.00003)})

test_that("analytic normal SW case matches a constant: nclusters", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = NA, 
                                        ntimes = 5, 
                                        d = 1.75, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = 0.5, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               300)})

test_that("analytic normal SW case matches a constant: ntimes", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = 300, 
                                        ntimes = NA, 
                                        d = 1.75, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = 0.5, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               5)})

test_that("analytic normal SW case matches a constant: d", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = 300, 
                                        ntimes = 5, 
                                        d = NA, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = 0.5, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               1.75)})

test_that("analytic normal SW case matches a constant: icc", {
  expect_equal(signif(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = 300, 
                                        ntimes = 5, 
                                        d = 1.75, 
                                        icc = NA, 
                                        rho_c = 0.5, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 3), 
               signif(clusterPower::createMissingVarianceParam(sigma_sq = c(100), 
                                                               sigma_b_sq = c(30)), 3))})

test_that("analytic normal SW case matches a constant: rho_c", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = 300, 
                                        ntimes = 5, 
                                        d = 1.75, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = NA, 
                                        rho_s = 0.25, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               0.4999818)})

test_that("analytic normal SW case matches a constant: rho_s", {
  expect_equal(as.numeric(cpa.sw.normal(nsubjects = 50.00003, 
                                        nclusters = 300, 
                                        ntimes = 5, 
                                        d = 1.75, 
                                        icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)), 
                                        rho_c = 0.4999818, 
                                        rho_s = NA, 
                                        vart = 2322.988, 
                                        power = 0.8100001)), 
               0.2502912)})