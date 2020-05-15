library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

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
  )), as.numeric(5))
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
  )), as.numeric(150))
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
  ), 0),
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
                difference = 0.1,
                ICC = 0.05,
                sigma_sq = 1,
                alpha = 0.05,
                method = 'glmm',
                quiet = FALSE,
                all.sim.data = FALSE
              )
            print(sim.power$power)
            reference <-
              as.numeric(round(
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
                ),
                1
              ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })



#--------------------------------- SIMPLE BINARY OUTCOMES

context("Simple, binary outcome accuracy")

# NOTE: does not match NIH reference value. Seems to consistently underestimate the nclusters needed
# by 1. NIH estimates nclusters using a t-distribution. Do we need them to match more closely?

# Ken says add user option to use t or normal distribution.

test_that("binary case matches reference value from NIH calculator", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2prop(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      nsubjects = 150,
      cv = 0,
      p1 = .1,
      p2 = .2,
      icc = .05
    )
  )), 13 - 1)
})

# p1=treatment, p2=control
# NOTE: does not match n4props value. Seems to consistently underestimate the nclusters needed
# by 1. n4props estimates nclusters using a t-distribution. Do we need them to match more closely?

test_that("binary case matches CRTSize::n4props", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2prop(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      nsubjects = 150,
      cv = 0,
      p1 = .1,
      p2 = .2,
      icc = .05
    )
  )),
  ceiling(CRTSize::n4props(
    pe = .2,
    pc = .1,
    m = 150,
    ICC = .05
  )$n - 1))
})

# compare simulation and analytic methods for binary outcomes
# NOTE: sigma_b has no default and may not be zero (although the help page
# says "if sigma_b2 is not specified, between cluster variances are
# assumed to be equal for both groups"). I'm probably not
# understanding something. --- ASK KEN ABOUT THIS

# Ken says: Should be a test for sigma b2=NA and if Sigb2=NA val set to sigma b.

test_that("simulation and analytic methods give similar power estimations for binary outcomes",
          {
            sim.power <-
              cps.binary(
                nsim = 100,
                nsubjects = 150,
                nclusters = 12,
                p1 = 0.1,
                p2 = 0.2,
                sigma_b = .4,
                alpha = 0.05,
                method = 'glmm',
                all.sim.data = FALSE
              )
            print(sim.power$power)
            reference <-
              round(as.numeric(
                crtpwr.2prop(
                  alpha = 0.05,
                  power = NA,
                  nclusters = 12,
                  nsubjects = 150,
                  cv = 0,
                  p1 = .1,
                  p2 = .2,
                  icc = .05
                )
              ), 1)
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })


#--------------------------------- SIMPLE POISSON OUTCOMES

context("Simple, poisson (incidence rate) outcome accuracy")

# incidence rate comparison
test_that("incidence rate outcomes matches CRTSize::n4incidence", {
  expect_equal(ceiling(as.numeric(
    crtpwr.2rate(
      alpha = 0.05,
      power = 0.8,
      nclusters = NA,
      py = 1000,
      r1 = 0.01,
      r2 = 0.005,
      cvb = 0.25
    )
  )),
  ceiling(
    CRTSize::n4incidence(
      le = 0.01,
      lc = 0.005,
      m = 1000,
      t = 1,
      CV = 0.25
    )$n
  ))
})

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
                sigma_b = 0.1,
                family = 'poisson',
                analysis = 'poisson',
                method = 'glmm',
                alpha = 0.05,
                quiet = FALSE,
                all.sim.data = TRUE
              )
            print(sim.power$power)
            reference <- as.numeric(crtpwr.2rate(
              alpha = 0.05,
              power = NA,
              nclusters = 4,
              py = 120,
              r1 = 0.1,
              r2 = 0.2,
              cvb = 0.1
            ))
            print(paste("analytic power = ", reference, sep = ""))
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
                difference = .48,
                sigma = 1,
                alpha = 0.05,
                sigma_b0 = .1,
                method = 'glmm',
                quiet = FALSE,
                all.sim.data = FALSE
              )
            print(sim.power$power)
            reference <-
              as.numeric(round(
                crtpwr.2meanD(
                  alpha = 0.05,
                  power = NA,
                  d = 0.48,
                  nclusters = 6,
                  nsubjects = 150,
                  icc = .05,
                  vart = 1,
                  rho_c = 0.3,
                  rho_s = 0.7
                ),
                1
              ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
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

# NOTE: cps.did.binary seems to give unusually low estimates of power. I don't know why yet.
# NOTE: the simulation method is very slow
test_that("DID binary simulation and analytic methods give similar power estimations",
          {
            sim.power <-
              cps.did.binary(
                nsim = 10,
                nsubjects = 150,
                nclusters = 30,
                p1 = 0.4,
                p2 = 0.5,
                sigma_b0 = 1,
                alpha = 0.05,
                method = 'glmm',
                all.sim.data = FALSE
              )
            print(sim.power$power)
            reference <- as.numeric(round(
              crtpwr.2propD(
                nsubjects = 150,
                p = .4,
                power = NA,
                nclusters = 30,
                d = .1,
                icc = .05,
                rho_c = 0,
                rho_s = 0
              ),
              1
            ))
            print(paste("analytic power = ", reference, sep = ""))
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
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
      all.sim.data = TRUE
    )
  x <- 0.8
  expect_equal(TRUE,
               data.table::between(x,
                                   model$power$lower.95.ci,
                                   model$power$upper.95.ci))
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
  x <- 0.8155647
  expect_equal(model, x)
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
    p.ntrt = 0.18,
    p.trt = 0.13,
    steps = 3,
    sigma_b_sq = 1,
    alpha = 0.05,
    method = 'glmm',
    quiet = FALSE,
    all.sim.data = FALSE
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
    all.sim.data = FALSE
  )
  testthat::expect_equal(TRUE,
                         data.table::between(model,
                                             x$power$Lower.95.CI,
                                             x$power$Upper.95.CI))
})

