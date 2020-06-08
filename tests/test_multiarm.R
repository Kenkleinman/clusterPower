library("devtools")
install_github("Kenkleinman/clusterPower@lexi")

library(testthat)
library(CRTSize)

library("clusterPower")

test_check("clusterPower")



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
                nsim = 1000,
                nsubjects = 50,
                narms = 2,
                nclusters = 4,
                probs = c(0.4, 0.1),
                sigma_b_sq = 1,
                alpha = 0.05,
                quiet = FALSE,
                method = 'glmm',
                all.sim.data = FALSE,
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
                sigma_b = 1,
                alpha = 0.05,
                method = 'glmm',
                all.sim.data = FALSE,
                seed = 123
              )[[3]][1]
            ))
          })

# compare simulation and analytic methods for continuous outcomes -uses ICC, asked John in include the definition in the documentation
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
                all.sim.data = FALSE
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
                all.sim.data = FALSE,
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
      all.sim.data = FALSE,
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
      all.sim.data = FALSE,
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
        all.sim.data = FALSE,
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
        all.sim.data = FALSE
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
            icc. <- runif(q, min = 0.1, max = 0.99)
            sig <- runif(q, min = 0.01, max = 2)
            for (i in 1:q) {
              multi.cps.normal <- cps.ma.normal(
                nsim = 100,
                narms = 2,
                nclusters = nc[i],
                nsubjects = ns[i],
                means = c(0, 1),
                ICC = icc.[i],
                sigma_sq = sig[i],
                alpha = 0.05,
                quiet = FALSE,
                method = 'glmm',
                all.sim.data = FALSE,
                poor.fit.override = TRUE,
                optmethod = "NLOPT_LN_NELDERMEAD",
                nofit = FALSE
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
                  all.sim.data = FALSE,
                  nofit = FALSE
                )
              rownames(multi.cps.normal$power) <- "Treatment.1"
              expect_equal(multi.cps.normal$power,
                           twoarm.mean$power)
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
      all.sim.data = FALSE,
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
      all.sim.data = TRUE
    )
    expect_equal(as.numeric(round(count.ma[[1]][, 1], 1)), as.numeric(round(count.sim$power[1], 1)))
    print(paste("Interation", i, "of 10."))
  } # end of loop
})
