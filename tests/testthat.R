#library(testthat)
#library(clusterPower)
#library(CRTSize)

#######################################
#### SIMPLE CRT DESIGNS TESTING #######
#######################################

#--------------------------------- SIMPLE CONTINUOUS OUTCOMES

context("Simple, normal outcome accuracy")

# compare to a reference value from NIH calculator
test_that("continuous case matches reference", {
  expect_equal(ceiling(as.numeric(
    cpa.normal(
      alpha = 0.05,
      power = 0.8,
      d = 0.4804988,
      nsubjects = 150,
      ICC = .05,
      vart = 1,
      method = "weighted",
      tol = .Machine$double.eps ^
        0.25
    )
  ))[1], as.numeric(5))
})

# compare to a reference value from PASS11
test_that("continuous case matches reference", {
  expect_equal(ceiling(as.numeric(
    cpa.normal(
      alpha = 0.05,
      power = 0.8,
      d = 0.4804988,
      nsubjects = NA,
      nclusters = 5,
      ICC = .05,
      vart = 1,
      method = "weighted",
      tol = .Machine$double.eps ^
        0.25
    )
  ))[1], as.numeric(150))
})


# compare to a comparable function in CRTSize
test_that("continuous case matches CRTSize", {
  expect_equal(round(as.numeric(
    cpa.normal(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      nsubjects = 150,
      ICC = .05,
      vart = 1,
      method = "taylor",
      tol = .Machine$double.eps ^
        0.25,
      d = 0.4804988
    )
  ), 0)[1],
  round(
    CRTSize::n4means(
      delta = 0.4804988,
      sigma = 1,
      m = 150,
      ICC = 0.05
    )$n
  ), 0)
})

# compare simulation and analytic methods for continuous outcomes
test_that("simulation and analytic methods give similar power estimations",
          {
            sim.power <-
              cps.normal(
                nsim = 1000,
                nsubjects = 100,
                nclusters = 100,
                mu = 1,
                mu2 = 1.1,
                ICC = 0.05,
                sigma_sq = 1,
                alpha = 0.05,
                method = 'glmm',
                quiet = FALSE,
                allSimData = FALSE
              )
            print(sim.power$power)
            reference <-
              as.numeric(
                cpa.normal(
                  alpha = 0.05,
                  power = NA,
                  d = 0.1,
                  nclusters = 100,
                  nsubjects = 100,
                  ICC = .05,
                  sigma_sq = 1,
                  method = "weighted",
                  tol = .Machine$double.eps ^ 0.25
                )
                )[1]
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })



#--------------------------------- SIMPLE BINARY OUTCOMES

context("Simple, binary outcome accuracy")

test_that("binary case matches reference value from NIH calculator", {
  expect_equal(round(as.numeric(
    cpa.binary(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      nsubjects = 150,
      CV = 0,
      p1 = .1,
      p2 = .2,
      ICC = .05,
      tdist = TRUE
    )
  ), 0), as.numeric(12))
})

test_that("binary case matches reference value from PASS 11", {
  expect_equal(round(as.numeric(
    cpa.binary(
      alpha = 0.05,
      power = NA,
      nclusters = 13,
      nsubjects = 150,
      CV = 0,
      p1 = .1,
      p2 = .2,
      ICC = .05,
      tdist = FALSE
    )
  ), 3), 0.86)
})


test_that("binary case matches CRTSize::n4props", {
  expect_equal(round(as.numeric(
    cpa.binary(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      nsubjects = 150,
      CV = 0,
      p1 = .1,
      p2 = .2,
      ICC = .05,
      tdist = FALSE
    )
  ), 0),
  round(CRTSize::n4props(
    pe = .2,
    pc = .1,
    m = 150,
    ICC = .05
  )$n - 1, 0))
})

# compare simulation and analytic methods for binary outcomes

# Ken says: Should be a test for sigma b2=NA and if Sigb2=NA val set to sigma b.
test_that("if sigma_b_sq2=NA, set to sigma_b_sq", {
  sim.power <-
    cps.binary(
      nsim = 100,
      nsubjects = 150,
      nclusters = 12,
      p1 = 0.1,
      p2 = 0.2,
      sigma_b_sq = .4,
      alpha = 0.05,
      method = 'glmm',
      allSimData = FALSE
    )
  expect_identical(sim.power$variance.parms[1], sim.power$variance.parms[2])
})


test_that("simulation and analytic methods give similar power estimations for binary outcomes",
          {
            sim.power <-
              cps.binary(
                nsim = 100,
                nsubjects = 150,
                nclusters = 12,
                p1 = 0.1,
                p2 = 0.2,
                sigma_b_sq = .4,
                alpha = 0.05,
                method = 'glmm',
                allSimData = FALSE
              )
            reference <-
              as.numeric(
                cpa.binary(
                  alpha = 0.05,
                  power = NA,
                  nclusters = 12,
                  nsubjects = 150,
                  CV = 0,
                  p1 = .1,
                  p2 = .2,
                  ICC = .05
                )
              )
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })


#--------------------------------- SIMPLE POISSON OUTCOMES

context("Simple, poisson (incidence rate) outcome accuracy")

# incidence rate comparison
test_that("incidence rate outcomes matches CRTSize::n4incidence", {
  expect_equal(round(as.numeric(
    cpa.count(
      alpha = 0.05,
      power = 0.8,
      nsubjects = 1000,
      nclusters = NA,
      r1 = 0.01,
      r2 = 0.005,
      CVB = 0.25
    )
  ), 0),
  round(
    CRTSize::n4incidence(
      le = 0.01,
      lc = 0.005,
      m = 1000,
      t = 1,
      CV = 0.25
    )$n, 0)
  )
})


#test_that("incidence rate outcomes matches PASS 11", {
#  expect_equal(round(as.numeric(
#    cpa.count(
#      alpha = 0.05,
#      power = 0.8,
#      nsubjects = 1000,
#      nclusters = NA,
#      r1 = 0.01,
#      r2 = 0.005,
#      CVB = 0.25
#    )
#  ), 0),
#  
#  )
#})

# compare simulation and analytic methods for poisson outcomes
# test .7-.9 in increments of .25
# NOTE: c1 & c2 (treatment) are outcome counts while r1 (treatment) & r2 are outcome rates
test_that("simulation and analytic methods give similar power estimations for poisson outcomes",
          {
            sim.power <-
              cps.count(
                nsim = 100,
                nsubjects = 120,
                nclusters = 4,
                c1 = 20,
                c2 = 10,
                sigma_b_sq = 0.1,
                family = 'poisson',
                analysis = 'poisson',
                method = 'glmm',
                alpha = 0.05,
                quiet = FALSE,
                allSimData = FALSE
              )
            reference <- as.numeric(cpa.count(
              alpha = 0.05,
              power = NA,
              nclusters = 4,
              nsubjects = 120,
              r1 = 0.1,
              r2 = 0.2,
              CVB = 0.1
            ))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })


#######################################
###### DID CRT DESIGNS TESTING ########
#######################################

#--------------------------------- DID CONTINUOUS OUTCOMES

context("DID normal outcome accuracy")

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cross-sectional)", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2meanD(
      alpha = 0.05,
      power = 0.8,
      d = 0.48,
      nsubjects = 150,
      rho_c = 0,
      rho_s = 0,
      icc = .05,
      vart = 1
    )
  )),
  as.numeric(9))
})

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2meanD(
      alpha = 0.05,
      power = 0.8,
      d = 0.48,
      nsubjects = 150,
      rho_c = 0.3,
      rho_s = 0.7,
      icc = .05,
      vart = 1
    )
  )),
  as.numeric(7))
})


test_that("DID normal simulation and analytic methods give similar power estimations",
          {
            sim.power <-
              cps.did.normal(
                nsim = 100,
                nsubjects = 150,
                nclusters = 6,
                mu = 1,
                mu2 = 1.48,
                sigma_sq = 1,
                alpha = 0.05,
                sigma_b_sq0 = .1,
                method = 'glmm',
                quiet = FALSE,
                allSimData = FALSE
              )
            print(sim.power$power)
            reference <-
              as.numeric(
                cpa.did.normal(
                  alpha = 0.05,
                  power = NA,
                  d = 0.48,
                  nclusters = 6,
                  nsubjects = 150,
                  ICC = .05,
                  vart = 1,
                  rho_c = 0.3,
                  rho_s = 0.7
                ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(TRUE,
                         data.table::between(reference,
                                             sim.power$power[, 2],
                                             sim.power$power[, 3]))
          })


#--------------------------------- DID BINARY OUTCOMES

context("DID binary outcome accuracy")

# compare to a reference value from NIH calculator
test_that("DID binary case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2propD(
      nsubjects = 50,
      p = .5,
      d = .1,
      icc = .05,
      rho_c = .3,
      rho_s = .7
    )
  )),
  as.numeric(33))
})


test_that("DID binary simulation and analytic methods give similar power estimations",
          {
            sim.power <-
              cps.did.binary(
                nsim = 10,
                nsubjects = 150,
                nclusters = 30,
                p1 = 0.4,
                p2 = 0.5,
                sigma_b_sq0 = 1,
                alpha = 0.05,
                method = 'glmm',
                allSimData = FALSE
              )
            print(sim.power$power)
            reference <- as.numeric(round(
              cpa.did.binary(
                nsubjects = 150,
                p = .4,
                power = NA,
                nclusters = 30,
                d = .1,
                ICC = .05,
                rho_c = 0,
                rho_s = 0
              ),
              1
            ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(TRUE,
                         data.table::between(reference,
                                             sim.power$power[, 2],
                                             sim.power$power[, 3]))
          })


#--------------------------------- DID COUNT OUTCOMES

context("DID count outcome accuracy")

# compare to a reference value from NIH calculator
test_that("DID count case matches reference (cohort)", {
  model <-
    cps.did.count(
      nsim = 100,
      nsubjects = 9,
      nclusters = 7,
      c1 = 5,
      c2 = 3,
      sigma_b_sq0 = c(1, 0.5),
      sigma_b_sq1 = c(0.5, 0.8),
      family = 'poisson',
      analysis = 'poisson',
      method = 'glmm',
      alpha = 0.05,
      quiet = FALSE,
      allSimData = TRUE,
      seed = 123
    )
  x <- 0.8
  expect_equal(TRUE,
               data.table::between(x,
                                   model$power$Lower.95.CI,
                                   model$power$Upper.95.CI))
})



##------------------------------------------ SW binary outcome

context("SW binary outcome accuracy")

# compare to a reference value from Zhou's FORTRAN program
test_that("SW binary case matches reference", {
  model <- cpa.sw.binary(
    nclusters = 20,
    steps = 3,
    nsubjects = 90,
    d = -0.05,
    ICC = 0.01,
    beta = -0.05,
    mu = 0.18
  )
  x <- round(0.8155647, 7)
  expect_equal(round(model, 7), x)
})

# compare SW binary analytic to simulated method
test_that("Analytic SW binary case matches simulation results", {
  model <- cpa.sw.binary(
    nclusters = 21,
    steps = 3,
    nsubjects = 90,
    d = -0.05,
    ICC = 0.01,
    beta = -0.05,
    mu = 0.18
  )
  x <- cps.sw.binary(
    nsim = 1000,
    nsubjects = 90,
    nclusters = 21,
    p1 = 0.18,
    p2 = 0.13,
    steps = 3,
    sigma_b_sq = 1,
    alpha = 0.05,
    method = 'glmm',
    quiet = FALSE,
    allSimData = FALSE,
    seed = 123
  )
  expect_equal(TRUE,
               data.table::between(model,
                                   x$power$Lower.95.CI,
                                   x$power$Upper.95.CI))
})

##------------------------------------------ SW count outcome
#compare analytic SW count to a constant

testthat::context("SW count outcome accuracy")

# compare to a reference value calculated previously
testthat::test_that("SW count case matches reference", {
  model <-
    cpa.sw.count(
      lambda1 = 1.75,
      RR = 0.9,
      nclusters = 21,
      steps = 6,
      nsubjects = 30,
      ICC = 0.01
    )
  x <- 0.806856
  testthat::expect_equal(round(model, 6), x)
})

# compare to a HH.count
testthat::test_that("SW count case matches HH.count", {
  model <-
    cpa.sw.count(
      lambda1 = 1.75,
      RR = 0.9,
      nclusters = 21,
      steps = 6,
      nsubjects = 30,
      ICC = 0.01
    )
  x <-
    HH.count(
      lambda1 = 1.75,
      RR = 0.9,
      I = 21,
      J = 6,
      K = 30,
      rho = 0.01
    )$power
  testthat::expect_equal(model, x)
})

# compare to simulated method
testthat::test_that("SW count case matches HH.count", {
  model <-
    cpa.sw.count(
      lambda1 = 1.75,
      RR = 0.9,
      nclusters = 25,
      steps = 5,
      nsubjects = 30,
      ICC = 0.01
    )
  x <- cps.sw.count(
    nsim = 1000,
    nsubjects = 30,
    nclusters = 25,
    c.ntrt = 1.75,
    c.trt = 1.575,
    steps = 5,
    sigma_b_sq = 0.01678129,
    alpha = 0.05,
    family = 'poisson',
    analysis = 'poisson',
    method = 'glmm',
    quiet = FALSE,
    allSimData = FALSE
  )
  testthat::expect_equal(TRUE,
                         data.table::between(model,
                                             x$power$Lower.95.CI,
                                             x$power$Upper.95.CI))
})



#######################################
######  IRGTT METHODS TESTING   #######
#######################################

#--------------------------------- IRGTT NORMAL OUTCOME

context("IRGTT outcome accuracy")

# compare analytic to simulation methods (normal)
test_that("analytic normal irgtt case matches simulated method normal irgtt case",
          {
            sim_model <- cps.irgtt.normal(
              nsim = 1000,
              nsubjects = c(400, 10),
              nclusters = 30,
              mu = 0,
              mu2 = 0.5,
              sigma_sq = 4,
              sigma_sq2 = 1,
              sigma_b_sq2 = 0.5,
              alpha = 0.05,
              quiet = FALSE,
              allSimData = FALSE,
              seed = 123
            )
            analytic_model <- cpa.irgtt.normal(
              nclusters = 30,
              nsubjects = 10,
              d = 0.5,
              varu = 0.5,
              vare = 1,
              varr = 4,
              ncontrols = 400,
              power = NA
            )
            expect_equal(
              TRUE,
              data.table::between(
                analytic_model,
                sim_model$power$Lower.95.CI,
                sim_model$power$Upper.95.CI
              )
            )
          })


# compare simulation method (normal) to a constant
test_that("simulation normal irgtt case matches a constant", {
  expect_equal(signif(as.numeric(
    cps.irgtt.normal(
      nsim = 1000,
      nsubjects = c(550, 12),
      nclusters = 125,
      mu = 0,
      mu2 = 0.25,
      sigma_sq = 3,
      sigma_sq2 = 2.1,
      sigma_b_sq2 = 0.25,
      alpha = 0.05,
      quiet = FALSE,
      allSimData = FALSE,
      seed = 123
    )$power$Power
  ), 2),
  0.78)
})

# compare analytic to a constant (normal)
test_that("analytic normal irgtt case matches a constant", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.normal(
      nclusters = 125,
      nsubjects = 12,
      d = 0.25,
      varu = 0.25,
      vare = 2.1,
      varr = 3,
      ncontrols = 550,
      power = NA
    ),
    2
  )),
  0.76)
})

#--------------------------------- IRGTT BINARY OUTCOME

# compare analytic to simulation methods (binary)
test_that("analytic binary irgtt power is within simulated method binary irgtt CI",
          {
            sim_model <- cps.irgtt.binary(
              nsim = 1000,
              nsubjects = 30,
              nclusters = 10,
              p1 = 0.2,
              p2 = 0.5,
              sigma_b_sq2 = 1,
              alpha = 0.05,
              allSimData = FALSE,
              seed = 123
            )
            analytic_model <- cpa.irgtt.binary(
              nclusters = 10,
              nsubjects = 30,
              ncontrols = 30,
              icc = 0.25,
              p.e = 0.5,
              p.c = 0.2,
              power = NA
            )
            expect_equal(
              TRUE,
              data.table::between(
                analytic_model,
                sim_model$power$Lower.95.CI,
                sim_model$power$Upper.95.CI
              )
            )
          })

# compare simulation methods (binary) to a constant
test_that("simulated method binary irgtt case matches a constant", {
  expect_equal(signif(as.numeric(
    cps.irgtt.binary(
      nsim = 1000,
      nsubjects = c(300, 12),
      nclusters = 120,
      p1 = 0.3,
      p2 = 0.2,
      sigma_b_sq2 = 1.75,
      alpha = 0.05,
      allSimData = FALSE,
      seed = 123
    )$power$Power
  ), 2),
  0.81)
})


# compare analytic methods (binary) to a constant
test_that("analytic binary irgtt case matches a constant: power", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      icc = 0.4,
      p.e = 0.2,
      p.c = 0.3,
      power = NA
    ),
    2
  )),
  0.79)
})

test_that("analytic binary irgtt case matches a constant: icc calc", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      icc = NA,
      p.e = 0.2,
      p.c = 0.3,
      power = 0.79
    ),
    2
  )),
  0.41)
})

test_that("analytic binary irgtt case matches a constant: nsubjects", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = NA,
      ncontrols = 300,
      icc = 0.4054797,
      p.e = 0.2,
      p.c = 0.3,
      power = 0.79
    ),
    2
  )),
  12)
})

test_that("analytic binary irgtt case matches a constant: nclusters", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = NA,
      nsubjects = 12,
      ncontrols = 300,
      icc = 0.4054797,
      p.e = 0.2,
      p.c = 0.3,
      power = 0.79
    ),
    2
  )),
  120)
})

test_that("analytic binary irgtt case matches a constant: p.c", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      icc = 0.4054797,
      p.e = 0.2,
      p.c = NA,
      power = 0.79
    ),
    2
  )),
  0.3)
})

test_that("analytic binary irgtt case matches a constant: p.e", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      icc = 0.4054797,
      p.e = NA,
      p.c = 0.3,
      power = 0.79
    ),
    2
  )),
  0.2)
})

#--------------------------------- IRGTT COUNT OUTCOME

# compare simulation methods (count) to a constant
test_that("analytic irgtt power is within simulated method count irgtt CI", {
  irgtt.count.sim <-
    cps.irgtt.count(
      nsim = 100,
      nsubjects = c(500, 10),
      nclusters = 100,
      c1 = 85,
      c2 = 450,
      sigma_b_sq2 = 25,
      family = 'poisson',
      analysis = 'poisson',
      alpha = 0.05,
      quiet = FALSE,
      allSimData = FALSE,
      opt = "L-BFGS-B"
    )
  
  sim_model <- cps.irgtt.count(
    nsim = 100,
    nsubjects = 30,
    nclusters = 10,
    p1 = 0.2,
    p2 = 0.5,
    sigma_b_sq2 = 1,
    alpha = 0.05,
    allSimData = FALSE,
    seed = 123
  )
  analytic_model <- cpa.irgtt.binary(
    nclusters = 10,
    nsubjects = 30,
    ncontrols = 30,
    icc = 0.25,
    p.e = 0.5,
    p.c = 0.2,
    power = NA
  )
  expect_equal(
    TRUE,
    data.table::between(
      analytic_model,
      sim_model$power$Lower.95.CI,
      sim_model$power$Upper.95.CI
    )
  )
})


###############################################
######  STEPPED-WEDGE METHODS TESTING   #######
###############################################

#--------------------------------- SW NORMAL OUTCOME

context("SW outcome accuracy")

# compare analytic to simulation methods (normal)
test_that("analytic normal SW outcome is within estimated CI of simulated SW method",
          {
            sim_model <- cps.sw.normal(
              nsim = 100,
              nsubjects = 50,
              nclusters = 30,
              difference = 1.75,
              steps = 5,
              sigma_sq = 100,
              sigma_b_sq = 30,
              alpha = 0.05,
              method = 'glmm',
              quiet = FALSE,
              allSimData = FALSE
            )
            analytic_model <- cpa.sw.normal(
              nsubjects = 50,
              nclusters = 300,
              ntimes = 5,
              d = 1.75,
              icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
              rho_c = 0.5,
              rho_s = 0.25,
              vart = 2322.988,
              power = NA
            )
            expect_equal(
              TRUE,
              data.table::between(
                analytic_model,
                sim_model$power$Lower.95.CI,
                sim_model$power$Upper.95.CI
              )
            )
          })

# compare SW simulation method to a constant (normal)
test_that("simulated method normal SW case matches a constant", {
  expect_equal(signif(as.numeric(
    cps.sw.normal(
      nsim = 100,
      nsubjects = 50,
      nclusters = 30,
      difference = 1.75,
      steps = 5,
      sigma_sq = 100,
      sigma_b_sq = 30,
      alpha = 0.05,
      method = 'glmm',
      quiet = FALSE,
      allSimData = FALSE
    )$power$Power
  ), 1),
  0.8)
})

# compare analytic SW methods (normal) to a constant
test_that("analytic normal SW case matches a constant", {
  expect_equal(as.numeric(signif(
    cpa.sw.normal(
      nsubjects = 50,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = NA
    ),
    1
  )),
  0.8)
})

test_that("analytic normal SW case matches a constant: vart", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = NA,
      power = 0.8100001
    )
  ),
  2322.988)
})

test_that("analytic normal SW case matches a constant: nsubjects", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = NA,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  50.00003)
})

test_that("analytic normal SW case matches a constant: nclusters", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = NA,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  300)
})

test_that("analytic normal SW case matches a constant: ntimes", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = NA,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  5)
})

test_that("analytic normal SW case matches a constant: d", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = 5,
      d = NA,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  1.75)
})

test_that("analytic normal SW case matches a constant: icc", {
  expect_equal(signif(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = NA,
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ), 3),
  signif(
    clusterPower::createMissingVarianceParam(sigma_sq = c(100),
                                             sigma_b_sq = c(30)),
    3
  ))
})

test_that("analytic normal SW case matches a constant: rho_c", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = NA,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  0.4999818)
})

test_that("analytic normal SW case matches a constant: rho_s", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      icc = clusterPower::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.4999818,
      rho_s = NA,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  0.2502912)
})


test_that("analytic normal SW case matches a function from SWSamp pkg", {
  expect_equal(as.numeric(
    cpa.sw.normal(
      alpha = 0.05,
      power = NA,
      nclusters = 4,
      nsubjects = 20,
      ntimes = 4,
      d = 5,
      icc = 0.05,
      rho_c = 1,
      rho_s = 0,
      vart = 400
    )
  ),
  as.numeric(
    SWSamp::HH.normal(
      mu = 0,
      b.trt = 5,
      sigma = sqrt(400),
      I = 4 * 4,
      J = 4,
      K = 20,
      rho = 0.05,
      which.var = "total"
    )$power
  ))
})


#######################################
#### MULTI-ARM CRT DESIGNS TESTING #######
#######################################

#--------------------------------- MULTI-ARM BINARY OUTCOME

context("Multi-arm outcome accuracy")

# compare to a reference value from NIH calculator
test_that("binary multi-arm case matches 2-arm binary case (simulated method)",
          {
            expect_equal(as.numeric(
              cps.ma.binary(
                nsim = 100,
                nsubjects = 50,
                narms = 2,
                nclusters = 4,
                probs = c(0.4, 0.1),
                sigma_b_sq = 1,
                alpha = 0.05,
                quiet = FALSE,
                method = 'glmm',
                allSimData = FALSE,
                multi.p.method = "none",
                poor.fit.override = TRUE,
                seed = 123,
                cores = "all"
              )[[1]][1]
            ),
            as.numeric(
              cps.binary(
                nsim = 100,
                nsubjects = 50,
                nclusters = 4,
                p1 = 0.4,
                p2 = 0.1,
                sigma_b_sq = 1,
                alpha = 0.05,
                method = 'glmm',
                allSimData = FALSE,
                seed = 123
              )[[3]][1]
            ))
          })

# compare simulation and analytic methods for continuous outcomes -uses ICC, 
#asked John in include the definition in the documentation
test_that("simulation and analytic methods give similar power estimations",
          {
            sim.power <-
              cps.normal(
                nsim = 1000,
                nsubjects = 10,
                nclusters = 16,
                difference = 1,
                ICC = 0.1009569,
                sigma_sq = 5,
                alpha = 0.05,
                method = 'glmm',
                quiet = FALSE,
                allSimData = FALSE
              )
            reference <-
              as.numeric(round(
                cpa.normal(
                  alpha = 0.05,
                  power = NA,
                  d = 1,
                  nclusters = 16,
                  nsubjects = 10,
                  ICC = 0.1009569,
                  vart = 5,
                  method = "weighted",
                  tol = .Machine$double.eps ^
                    0.25
                ),
                1
              ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })

# compare to a comparable function in CRTSize -uses ICC
test_that("continuous case matches CRTSize", {
  nc <- sample.int(200, 10)
  ns <- sample.int(200, 10)
  icc. <- runif(10, min = 0.1, max = 0.99)
  for (i in 1:10) {
    expect_equal(round(as.numeric(
      cps.normal(
        alpha = 0.05,
        power = 0.8,
        nclusters = nc[i],
        nsubjects = ns[i],
        ICC = icc.[i],
        vart = 1,
        method = "weighted",
        tol = .Machine$double.eps ^
          0.25,
        d = 1
      )
    ), 0),
    round(n4means(
      delta = 1,
      sigma = sqrt(1),
      m = 10,
      ICC = icc.[i]
    )$n, 0))
    print(nc[i])
    print(ns[i])
    print(icc.[i])
    print(sig[i])
    print(paste("Interation", i, "of 10."))
  } # end of loop
})

#Normal simulation methods

test_that("continuous simulation method matches a reference (previous value)",
          {
            nsubjects.example <- list(
              c(20, 20, 20, 20, 20, 75, 20, 20, 20, 75),
              c(20, 20, 25, 25, 25, 25, 25, 25),
              c(40, 25, 40, 20, 20, 20, 20, 20)
            )
            means.example <- c(1, 1.75, 0)
            sigma_sq.example <- c(2, 1.2, 2)
            sigma_b_sq.example <- c(1.1, 1.15, 1.1)
            multi.cps.normal.unbal <-
              cps.ma.normal(
                nsim = 100,
                nsubjects = nsubjects.example,
                means = means.example,
                sigma_sq = sigma_sq.example,
                sigma_b_sq = sigma_b_sq.example,
                alpha = 0.05,
                quiet = FALSE,
                ICC = NULL,
                method = 'glmm',
                allSimData = FALSE,
                seed = 123,
                cores = "all",
                poor.fit.override = FALSE,
                opt = "nlminb"
              )
            prev.value <- t(data.frame(0.34, 0.45))
            prev.value <- data.frame(as.numeric(prev.value))
            rownames(prev.value) <- c("Treatment.2", "Treatment.3")
            colnames(prev.value) <- "Power"
            expect_equal(round(multi.cps.normal.unbal[['power']]['Power'], 2), prev.value)
          })


test_that("normal vs t-dist comparison", {
  q <- 10
  nc <- sample.int(200, q)
  ns <- sample.int(200, q)
  icc. <- runif(q, min = 0.1, max = 0.99)
  sig <- (runif(q, min = 0.1, max = 50)) ^ 2
  sigb <- createMissingVarianceParam(ICC = icc., sigma_b = sig)
  holder <- data.frame(nc, ns, icc., sig, sigb)
  same <- rep(NA, length = q)
  for (i in 1:q) {
    multi.cps.normal <- cps.ma.normal(
      nsim = 200,
      narms = 2,
      nclusters = nc[i],
      nsubjects = ns[i],
      means = c(0, 1),
      tdist = FALSE,
      ICC = icc.[i],
      sigma_sq = sig[i],
      alpha = 0.05,
      quiet = FALSE,
      method = 'glmm',
      allSimData = FALSE,
      poor.fit.override = TRUE,
      cores = "all"
    )
    multi.cps.tdist <- cps.ma.normal(
      nsim = 200,
      narms = 2,
      nclusters = nc[i],
      nsubjects = ns[i],
      means = c(0, 1),
      tdist = TRUE,
      ICC = icc.[i],
      sigma_sq = sig[i],
      alpha = 0.05,
      quiet = FALSE,
      method = 'glmm',
      allSimData = FALSE,
      poor.fit.override = TRUE,
      cores = "all"
    )
    if (round(multi.cps.normal[['power']]['Power'], 2) == round(multi.cps.tdist[['power']]['Power'], 2)) {
      same[i] <- 1
    }  else {
      same[i] <- 0
    }
  } # end of loop
  expect_equal(same, rep(1, times = q))
})

test_that("continuous simulation method matches the analytic method", {
  q <- 10
  nc <- sample.int(200, q)
  ns <- sample.int(200, q)
  icc. <- runif(q, min = 0.01, max = 0.99)
  sig <- runif(q, min = 0.01, max = 3)
  sigb <- createMissingVarianceParam(ICC = icc., sigma_b = sig)
  holder <- data.frame(nc, ns, icc., sig, sigb)
  for (i in 1:q) {
    multi.cps.normal <-
      cps.ma.normal(
        nsim = 100,
        narms = 2,
        nclusters = nc[i],
        nsubjects = ns[i],
        means = c(0, 1),
        tdist = FALSE,
        ICC = icc.[i],
        sigma_sq = sig[i],
        alpha = 0.05,
        quiet = FALSE,
        method = 'glmm',
        allSimData = FALSE,
        poor.fit.override = TRUE,
        low.power.override = TRUE,
        cores = NULL,
        optmethod = "NLOPT_LN_NELDERMEAD"
      )
    twoarm.mean <-
      cps.normal(
        nsim = 100,
        nsubjects = ns[i],
        nclusters = nc[i],
        mu = 0,
        mu2 = 1,
        ICC = icc.[i],
        sigma_sq = sig[i],
        alpha = 0.05,
        method = 'glmm',
        quiet = FALSE,
        allSimData = FALSE
      )
    analytic.mean <- cpa.normal(
      alpha = 0.05,
      power = NA,
      nclusters = nc[i],
      nsubjects = ns[i],
      ICC = icc.[i],
      sigma_sq = sig[i],
      method = "weighted",
      tol = .Machine$double.eps ^ 0.25,
      d = 1
    )
    print(paste("Interation", i, "of 10."))
    expect_equal(
      data.table::between(
        as.numeric(analytic.mean),
        twoarm.mean$power$Lower.95.CI,
        twoarm.mean$power$Upper.95.CI
      ),
      TRUE
    )
    expect_equal(
      data.table::between(
        as.numeric(analytic.mean),
        multi.cps.normal$power$Lower.95.CI,
        multi.cps.normal$power$Upper.95.CI
      ),
      TRUE
    )
    expect_equal(
      data.table::between(
        multi.cps.normal$power$Power,
        twoarm.mean$power$Lower.95.CI,
        twoarm.mean$power$Upper.95.CI
      ),
      TRUE
    )
    expect_equal(
      data.table::between(
        twoarm.mean$power$Power,
        multi.cps.normal$power$Lower.95.CI,
        multi.cps.normal$power$Upper.95.CI
      ),
      TRUE
    )
  } # end of loop
})


test_that("continuous simulation method matches the 2-arm simulation method",
          {
            q <- 1
            nc <- sample.int(200, q)
            ns <- sample.int(200, q)
            prob1 <- runif(q, min = 0, max = 1)
            prob2 <- runif(q, min = 0, max = 1)
            icc. <- runif(q, min = 0.1, max = 0.99)
            sig <- runif(q, min = 0.01, max = 2)
            for (i in 1:q) {
              multi.cps.bin <- cps.ma.binary(
                nsim = 100, nsubjects = ns[i], 
                narms = 2,
                nclusters = nc[i],
                probs = c(prob1[i], prob2[i]),
                sigma_b_sq = sig[i], alpha = 0.05,
                quiet = FALSE, method = 'glmm', 
                allSimData = FALSE, 
                multi.p.method = "none",
                seed = 123, cores = "all", 
                poor.fit.override = TRUE  
              )
              twoarm.mean <-
                cps.normal(
                  nsim = 100,
                  nsubjects = ns[i],
                  nclusters = nc[i],
                  mu = 0,
                  mu2 = 1,
                  ICC = icc.[i],
                  sigma_sq = sig[i],
                  alpha = 0.05,
                  method = 'glmm',
                  quiet = FALSE,
                  allSimData = FALSE,
                  nofit = FALSE
                )
            expect_equal(
              data.table::between(
                twoarm.mean$power$Power,
                multi.cps.normal$power$Lower.95.CI,
                multi.cps.normal$power$Upper.95.CI
              ),
              TRUE
            )
            print(paste("Interation", i, "of 10."))
          } # end of loop
          })


###### Multi-arm Count outcome testing

test_that("count simulation method matches the 2-arm simulation method", {
  q <- 1
  nc <- sample.int(200, q)
  ns <- sample.int(200, q)
  sig <- runif(q, min = 0.01, max = 1)
  for (i in 1:q) {
    count.ma <- cps.ma.count(
      nsim = 100,
      nsubjects = ns[i],
      narms = 2,
      nclusters = nc[i],
      counts = c(30, 70),
      sigma_b_sq = sig[i],
      alpha = 0.05,
      quiet = TRUE,
      method = 'glmm',
      allSimData = FALSE,
      multi.p.method = "none",
      poor.fit.override = TRUE
    )
    count.sim = cps.count(
      nsim = 100,
      nsubjects = ns[i],
      nclusters = nc[i],
      c1 = 30,
      c2 = 70,
      sigma_b = sig[i],
      family = 'poisson',
      analysis = 'poisson',
      method = 'glmm',
      alpha = 0.05,
      quiet = FALSE,
      allSimData = TRUE
    )
    expect_equal(as.numeric(round(count.ma[[1]][, 1], 1)), as.numeric(round(count.sim$power[1], 1)))
    print(paste("Interation", i, "of 10."))
  } # end of loop
})



############################################
######  TESTING FOR SMALL FUNCTIONS  #######
############################################

context("confintCalc")

test_that("confintCalc matches a reference", {
  some_nums <- confintCalc(
    nsim = 1000,
    alpha = 0.05,
    p.val = 0.01,
    names.power = c("bestCI")
  )
  expect_equal(3, length(some_nums))
  expect_equal(0.001, some_nums$Power)
  expect_equal(2.531749e-05, some_nums$Lower.95.CI)
  expect_equal(0.005558924, some_nums$Upper.95.CI)
})

context("createMissingVarianceParam")

test_that("createMissingVarianceParam matches a reference", {
  ICC <-
    createMissingVarianceParam(sigma_sq = c(1, 1, 0.9),
                               sigma_b_sq = c(0.1, 0.15, 0.1))
  expect_equal(3, length(ICC))
  expect_equal(0.09090909, ICC[1])
  expect_equal(0.13043478, ICC[2])
  expect_equal(0.10000000, ICC[3])
  sig_b <-
    createMissingVarianceParam(sigma_sq = c(1, 1, 0.9), ICC = ICC)
  expect_equal(3, length(sig_b))
  expect_equal(0.10, sig_b[1])
  expect_equal(0.15, sig_b[2])
  expect_equal(0.10, sig_b[3])
})

context("is.wholenumber")

test_that("is.wholenumber matches a reference", {
  expect_equal(FALSE, is.wholenumber(3.5))
  expect_equal(TRUE, is.wholenumber(5001))
})


context("prop_H0_rejection")

test_that("prop_H0_rejection matches a reference", {
  prop <- prop_H0_rejection(alpha = 0.05,
                            nsim = 1000,
                            LRT.holder.abbrev = 804)
  expect_equal(3, length(prop))
  expect_equal(0.804, as.numeric(prop[1]))
  expect_equal(0.828, as.numeric(prop[3]))
  expect_equal("Power", names(prop[1]))
})


context("type1ErrTest")

test_that("type1ErrTest matches a reference", {
  warn <- type1ErrTest(
    sigma_sq_ = c(0.1, 4),
    sigma_b_sq_ = c(0.1, 0.15),
    nsubjects_ = list(rep(4, 10),
                      rep(4, 10))
  )
  expect_equal(310, nchar(warn))
  nowarn <- type1ErrTest(
    sigma_sq_ = c(1, 1),
    sigma_b_sq_ = c(0.5, 0.1),
    nsubjects_ = list(rep(20, 20),
                      rep(20, 20))
  )
  expect_equal(TRUE, is.null(nowarn))
})

#optimizerSearch

context("optimizerSearch")

test_that("optimizerSearch matches a reference", {
  library(lattice)
  gm1 <-
    glmer(
      cbind(incidence, size - incidence) ~ period + (1 | herd),
      data = cbpp,
      family = binomial
    )
  lm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  lmopt <- optimizerSearch(lm1)
  gmopt <- optimizerSearch(gm1)
  expect_equal("bobyqa", gmopt)
  expect_equal("bobyqa", lmopt)
})

testthat::context(
  "createMissingVarianceParam: calculate ICC, sigma,
                  or sigma_b based on user providing 2 of those values"
)

# compare to a reference value
testthat::test_that("createMissingVarianceParam returns ICC", {
  testthat::expect_equal(
    createMissingVarianceParam(
      sigma = c(1, 1, 0.9),
      sigma_b = c(0.1, 0.15, 0.1)
    ),
    c(0.09090909, 0.13043478, 0.10000000)
  )
})

testthat::test_that("createMissingVarianceParam fails when fewer than 2 params provided",
                    {
                      testthat::show_failure(createMissingVarianceParam(sigma = c(1, 1, 0.9)))
                    })
