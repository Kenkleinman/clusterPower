#library(testthat)
#library(clusterPower)
#library(CRTSize)

#######################################
#### SIMPLE CRT DESIGNS TESTING #######
#######################################

#--------------------------------- SIMPLE CONTINUOUS OUTCOMES

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
  ), 0),
  round(
    CRTSize::n4means(
      delta = 0.4804988,
      sigma = 1,
      m = 150,
      ICC = 0.05
    )$n, 0)
  )
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
            expect_equal(sim.power$power[, 2] <= reference &
                           sim.power$power[, 3] >= reference,
                         TRUE)
          })



#--------------------------------- SIMPLE BINARY OUTCOMES

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

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cross-sectional)", {
  expect_equal(ceiling(as.numeric(
    cpa.did.normal(
      alpha = 0.05,
      power = 0.8,
      d = 0.48,
      nsubjects = 150,
      rho_c = 0,
      rho_s = 0,
      ICC = .05,
      vart = 1
    )
  )),
  as.numeric(9))
})

# compare to a reference value from NIH calculator
test_that("DID continuous case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(
    cpa.did.normal(
      alpha = 0.05,
      power = 0.8,
      d = 0.48,
      nsubjects = 150,
      rho_c = 0.3,
      rho_s = 0.7,
      ICC = .05,
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
                delta_mu2 = 0.48,
                sigma_sq = 0.8,
                alpha = 0.05,
                sigma_b_sq0 = 0.04210,
                method = 'glmm',
                quiet = FALSE,
                allSimData = FALSE,
                seed = 123
              )
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
            expect_equal(TRUE,
                         data.table::between(reference,
                                             sim.power$power[, 2],
                                             sim.power$power[, 3]))
          })


#--------------------------------- DID BINARY OUTCOMES

# compare to a reference value from NIH calculator
test_that("DID binary case matches reference (cohort)", {
  expect_equal(ceiling(as.numeric(
    cpa.did.binary(
      nsubjects = 50,
      p = .5,
      d = .1,
      ICC = .05,
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
                p1t0 = 0,
                p2t1 = 0.1,
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

# compare to a reference value from NIH calculator
test_that("DID count case matches reference (cohort)", {
  model <-
    cps.did.count(
      nsim = 100,
      nsubjects = 9,
      nclusters = 7,
      c1t0 = 5,
      c2t1 = 3,
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

# compare to a reference value from Zhou's FORTRAN program
test_that("SW binary case matches reference", {
  model <- cpa.sw.binary(nclusters = 9,
     steps = 3,
     nsubjects = 20,
     timeEffect = 0,
     ICC = 0.05,
     p1 = 0.31,
     p0 = 0.2)
  dataset <- matrix(c(rep(c(0,1,1,1),3),rep(c(0,0,1,1),3), rep(c(0,0,0,1), 3)),9,4,byrow=TRUE)
  x <- swdpwr::swdpower(K = 20, 
                        design = dataset, 
                        family = "binomial", 
                        model = "conditional", 
                        link = "identity", 
                        type = "cross-sectional", 
                        meanresponse_start = 0.2, 
                        meanresponse_end0 = 0.2, 
                        meanresponse_end1 = 0.31, 
                        typeIerror = 0.05, 
                        alpha0 = 0.05,  
                        alpha1 = 0.05)$Summary[["Power",1]]
  expect_equal(as.numeric(round(model, 3)), as.numeric(x))
})

# compare SW binary analytic to simulated method
test_that("Analytic SW binary case matches simulation results", {
  model <- cpa.sw.binary(
    nclusters = 12,
    steps = 3,
    nsubjects = 10,
    timeEffect = 0,
    ICC = 0.01,
    p1 = 0.2,
    p0 = 0.1
  )
  x <- cps.sw.binary(nsim = 10, 
                     nsubjects = 10, 
                     nclusters = 12, 
                     p0 = 0.1, 
                     p1 = 0.2, 
                     steps = 3, 
                     sigma_b_sq = 1, 
                     alpha = 0.05, 
                     method = 'glmm', 
                     quiet = FALSE, 
                     allSimData = FALSE, 
                     seed = 123)
   expect_equal(TRUE,
               data.table::between(model,
                                   x$power$Lower.95.CI,
                                   x$power$Upper.95.CI))
})

##------------------------------------------ SW count outcome
#compare analytic SW count to a constant

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
    SWSamp::HH.count(
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
testthat::test_that("SW analytic count case matched simulated outcome", {
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
    c0 = 1.75,
    c1 = 1.575,
    steps = 5,
    sigma_b_sq = 0.01678129,
    alpha = 0.05,
    family = 'poisson',
    analysis = 'poisson',
    method = 'glmm',
    quiet = FALSE,
    allSimData = FALSE,
    seed = 123
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
  ), 3),
  0.736)
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
              sigma_b_sq2 = 0.1,
              alpha = 0.05,
              allSimData = FALSE,
              seed = 123
            )
            analytic_model <- cpa.irgtt.binary(
              nclusters = 10,
              nsubjects = 30,
              ncontrols = 30,
              ICC = 0.2,
              p2 = 0.5,
              p1 = 0.2,
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
  ), 3),
  0.675)
})


# compare analytic methods (binary) to a constant
test_that("analytic binary irgtt case matches a constant: power", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      ICC = 0.4,
      p2 = 0.2,
      p1 = 0.3,
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
      ICC = NA,
      p2 = 0.2,
      p1 = 0.3,
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
      ICC = 0.4054797,
      p2 = 0.2,
      p1 = 0.3,
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
      ICC = 0.4054797,
      p2 = 0.2,
      p1 = 0.3,
      power = 0.79
    ),
    2
  )),
  120)
})

test_that("analytic binary irgtt case matches a constant: p1", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      ICC = 0.4054797,
      p2 = 0.2,
      p1 = NA,
      power = 0.79
    ),
    2
  )),
  0.3)
})

test_that("analytic binary irgtt case matches a constant: p2", {
  expect_equal(as.numeric(signif(
    cpa.irgtt.binary(
      nclusters = 120,
      nsubjects = 12,
      ncontrols = 300,
      ICC = 0.4054797,
      p2 = NA,
      p1 = 0.3,
      power = 0.79
    ),
    2
  )),
  0.2)
})

#--------------------------------- IRGTT COUNT OUTCOME

# compare simulation methods (count) to a constant
test_that("previous value (constant) is within simulated method count irgtt CI", {
  sim_model <- cps.irgtt.count(
    nsim = 100,
    nsubjects = 30,
    nclusters = 10,
    c1 = 1,
    c2 = 3,
    sigma_b_sq2 = 1,
    alpha = 0.05,
    allSimData = FALSE,
    seed = 123
  )
  expect_equal(
    TRUE,
    data.table::between(
      0.84,
      sim_model$power$Lower.95.CI,
      sim_model$power$Upper.95.CI
    )
  )
})


###############################################
######  STEPPED-WEDGE METHODS TESTING   #######
###############################################

#--------------------------------- SW NORMAL OUTCOME

# compare analytic to simulation methods (normal)
test_that("analytic normal SW outcome is within estimated CI of simulated SW method",
          {
            sim_model <- cps.sw.normal(
              nsim = 100,
              nsubjects = 50,
              nclusters = 30,
              mu1 = 1.75,
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
              ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      mu1 = 1.75,
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ),
  1.75)
})

test_that("analytic normal SW case matches a constant: ICC", {
  expect_equal(signif(as.numeric(
    cpa.sw.normal(
      nsubjects = 50.00003,
      nclusters = 300,
      ntimes = 5,
      d = 1.75,
      ICC = NA,
      rho_c = 0.5,
      rho_s = 0.25,
      vart = 2322.988,
      power = 0.8100001
    )
  ), 3),
  signif(
    clusterPower:::createMissingVarianceParam(sigma_sq = c(100),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = clusterPower:::createMissingVarianceParam(sigma_sq = c(100), sigma_b_sq = c(30)),
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
      ICC = 0.05,
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

#--------------------------------- MULTI-ARM NORMAL OUTCOME

# analytic

test_that("normal simulation methods match a reference I",
          {
          x <- cpa.ma.normal(narms=3,nclusters=5,vara=36,varc=9,vare=64)
            expect_equal(as.numeric(round(x, 0)), 21)
          })

# compare simulation and analytic methods for continuous outcomes -uses ICC, 
#asked John in include the definition in the documentation
test_that("normal simulation methods match a reference I",
          {
            nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
            means.example <- c(22, 21, 21.5)
            sigma_sq.example <- c(1, 1, 0.9)
            sigma_b_sq.example <- c(0.1, 0.15, 0.1)
            
            multi.cps.normal.models <- cps.ma.normal(nsim = 100,
                                                     narms = 3,
                                                     nsubjects = nsubjects.example,
                                                     means = means.example,
                                                     sigma_sq = sigma_sq.example,
                                                     sigma_b_sq = sigma_b_sq.example,
                                                     alpha = 0.05,
                                                     quiet = FALSE, method = 'glmm',
                                                     seed = 123, cores = "all",
                                                     lowPowerOverride = FALSE,
                                                     poorFitOverride = FALSE,
                                                     optmethod = "nlm")
            expect_equal(multi.cps.normal.models$power$Power, 0.84)
          })

test_that("normal simulation methods match a reference II",
          {
            multi.cps.normal <- cps.ma.normal(nsim = 100, narms = 3,
                                              nclusters = c(10,11,10), nsubjects = 25,
                                              means = c(1, 0.25, 1.75),
                                              sigma_sq = c(1.2, 1, 1.9),
                                              sigma_b_sq = c(0.5, 1, 0.75),
                                              quiet = FALSE, ICC=NULL, method = 'glmm',
                                              allSimData = FALSE, seed = 123,
                                              poorFitOverride = TRUE,
                                              cores = NULL,
                                              optmethod = "nlminb")
            expect_equal(multi.cps.normal$power$Power, 0.88)
          })

test_that("normal simulation methods match a reference III",
          {
 multi.cps.normal.simple <- cps.ma.normal(nsim = 1000, narms = 3,
                                   nclusters = 5, nsubjects = 10,
                                   means = c(22.0, 22.5, 22.9),
                                   sigma_sq = 0.2,
                                   sigma_b_sq = 0.2, alpha = 0.05,
                                   quiet = FALSE, ICC=NULL, method = 'glmm',
                                   allSimData = FALSE, seed = 123,
                                   poorFitOverride = TRUE, cores="all",
                                   optmethod = "NLOPT_LN_NELDERMEAD")
 expect_equal(multi.cps.normal.simple$power$Power, 0.781)
          })

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
                narms = 3,
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
                poorFitOverride = FALSE,
                opt = "nlminb"
              )
            prev.value <- t(data.frame(0.34, 0.45))
            prev.value <- data.frame(as.numeric(prev.value))
            rownames(prev.value) <- c("Arm.2", "Arm.3")
            colnames(prev.value) <- "Power"
            expect_equal(round(multi.cps.normal.unbal[['overall.power2']]['Power'], 2), 
                         prev.value)
          })

###### Multi-arm Binary outcome testing

test_that("continuous simulation method matches a reference (previous value) I",
          {
 bin.ma.rct.unbal <- cps.ma.binary(nsim = 12,
                             nsubjects = list(rep(20, times=15),
                             rep(15, times=15),
                             rep(17, times=15)),
                             narms = 3,
                             nclusters = 15,
                             probs = c(0.35, 0.43, 0.50),
                             sigma_b_sq = c(0.1, 0.1, 0.1),
                             alpha = 0.05, allSimData = TRUE,
                             seed = 123, cores="all")
 expect_equal(round(bin.ma.rct.unbal$power$Power, 6), 0.833333)
          })

test_that("continuous simulation method matches a reference (previous value) II",
          {
 bin.ma.rct.bal <- cps.ma.binary(nsim = 100, nsubjects = 50, narms=3,
                             nclusters = 8,
                             probs = c(0.35, 0.4, 0.5),
                             sigma_b_sq = 0.1, alpha = 0.05,
                             quiet = FALSE, method = 'glmm',
                             allSimData = FALSE,
                             multi_p_method="none",
                             seed = 123, cores="all",
                             poorFitOverride = FALSE)
 expect_equal(round(bin.ma.rct.bal$power$Power, 6), 0.808081)
          })

###### Multi-arm Count outcome testing

test_that("count simulation method matches a reference I", {
   nsubjects.example <- list(c(150, 200, 50, 100), c(50, 150, 210, 100),
                          c(70, 200, 150, 50, 100))
   counts.example <- c(10, 55, 65)
   sigma_b_sq.example <- c(1, 1, 2)
   count.ma.rct.unbal <- cps.ma.count(nsim = 100,
                               narms = 3,
                               nsubjects = nsubjects.example,
                               counts = counts.example,
                               sigma_b_sq = sigma_b_sq.example,
                               alpha = 0.05, seed = 123)
   expect_equal(round(count.ma.rct.unbal$power$Power, 6), 0.808081)
})
             
test_that("count simulation method matches a reference II", {                  
 count.ma.rct.bal <- cps.ma.count(nsim = 10, nsubjects = 100, narms = 4,
                               nclusters = 25, counts = c(30, 35, 70, 40),
                               sigma_b_sq = c(1, 1.2, 1, 0.9), seed = 123)
 expect_equal(round(count.ma.rct.bal$power$Power, 6), 0.9)
})



############################################
######  TESTING FOR SMALL FUNCTIONS  #######
############################################

test_that("confintCalc matches a reference", {
  some_nums <- clusterPower:::confintCalc(
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


test_that("createMissingVarianceParam matches a reference", {
  ICC <-
    clusterPower:::createMissingVarianceParam(sigma_sq = c(1, 1, 0.9),
                               sigma_b_sq = c(0.1, 0.15, 0.1))
  expect_equal(3, length(ICC))
  expect_equal(0.09090909, ICC[1])
  expect_equal(0.13043478, ICC[2])
  expect_equal(0.10000000, ICC[3])
  sig_b <-
    clusterPower:::createMissingVarianceParam(sigma_sq = c(1, 1, 0.9), ICC = ICC)
  expect_equal(3, length(sig_b))
  expect_equal(0.10, sig_b[1])
  expect_equal(0.15, sig_b[2])
  expect_equal(0.10, sig_b[3])
})


test_that("is.wholenumber matches a reference", {
  expect_equal(FALSE, clusterPower:::is.wholenumber(3.5))
  expect_equal(TRUE, clusterPower:::is.wholenumber(5001))
})



test_that("prop_H0_rejection matches a reference", {
  prop <- clusterPower:::prop_H0_rejection(alpha = 0.05,
                            nsim = 1000,
                            LRT.holder.abbrev = 804)
  expect_equal(3, length(prop))
  expect_equal(0.804, as.numeric(prop[1]))
  expect_equal(0.828, as.numeric(prop[3]))
  expect_equal("Power", names(prop[1]))
})



test_that("type1ErrTest matches a reference", {
  warn <- clusterPower:::type1ErrTest(
    sigma_sq_ = c(0.1, 4),
    sigma_b_sq_ = c(0.1, 0.15),
    nsubjects_ = list(rep(4, 10),
                      rep(4, 10))
  )
  expect_equal(310, nchar(warn))
  nowarn <- clusterPower:::type1ErrTest(
    sigma_sq_ = c(1, 1),
    sigma_b_sq_ = c(0.5, 0.1),
    nsubjects_ = list(rep(20, 20),
                      rep(20, 20))
  )
  expect_equal(TRUE, is.null(nowarn))
})

#optimizerSearch


test_that("optimizerSearch matches a reference", {
  gm1 <-
    lme4::glmer(
      cbind(incidence, size - incidence) ~ period + (1 | herd),
      data = cbpp,
      family = binomial
    )
  lm1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  lmopt <- clusterPower:::optimizerSearch(lm1)
  gmopt <- clusterPower:::optimizerSearch(gm1)
  expect_equal("bobyqa", gmopt)
  expect_equal("bobyqa", lmopt)
})


# compare to a reference value
testthat::test_that("createMissingVarianceParam returns ICC", {
  testthat::expect_equal(
    clusterPower:::createMissingVarianceParam(
      sigma = c(1, 1, 0.9),
      sigma_b = c(0.1, 0.15, 0.1)
    ),
    c(0.09090909, 0.13043478, 0.10000000)
  )
})

testthat::test_that("createMissingVarianceParam fails when fewer than 2 params provided",
                    {
                      testthat::show_failure(clusterPower:::createMissingVarianceParam(sigma = c(1, 1, 0.9)))
                    })
