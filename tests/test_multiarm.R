library(testthat)
library(clusterPower)
library(CRTSize)

test_check("clusterPower")

#######################################
#### MULTI-ARM CRT DESIGNS TESTING #######
#######################################

#--------------------------------- MULTI-ARM CONTINUOUS OUTCOMES

context("Multi-arm outcome accuracy")

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

# compare simulation and analytic methods for continuous outcomes -uses ICC, asked John in include the definition in the documentation
test_that("simulation and analytic methods give similar power estimations", {
  sim.power <- cps.normal(nsim = 100, nsubjects = 10, nclusters = 16, difference = 1,
                          ICC = 0.1009569, sigma = 5, alpha = 0.05, method = 'glmm', 
                          quiet = FALSE, all.sim.data = FALSE)
  reference <- as.numeric(round(crtpwr.2mean(alpha = 0.05, power = NA, d = 1, 
                                             nclusters = 16, nsubjects = 10,
                                             icc = 0.1009569, vart = 5, method = "weighted", 
                                             tol = .Machine$double.eps^0.25), 1))
  print(paste("analytic power = ", reference, sep=""))
  expect_equal(sim.power$power[,2] <= reference & sim.power$power[,3] >= reference, TRUE)
})

# FIXME fix the grid values and insert
# compare to a comparable function in CRTSize -uses ICC
test_that("continuous case matches CRTSize", {
  nc <- sample.int(200, 10)
  ns <- sample.int(200, 10)
  icc. <- runif(10, min=0.1, max=0.99)
  sig <- runif(10, min=0.1, max=50)
  for (i in 1:10){
  expect_equal(round(as.numeric(crtpwr.2mean(alpha = 0.05, power = 0.8, 
                                             nclusters = NA, nsubjects = 10,
                                             icc = 0.1009569, vart = 5, method = "weighted", 
                                             tol = .Machine$double.eps^0.25, d=1)), 0), 
               round(n4means(delta=1, sigma=sqrt(5), m=10, ICC=0.1009569)$n,0))
    print(nc[i])
    print(ns[i])
    print(icc.[i])
    print(sig[i])
    print(paste("Interation", i, "of 10."))
  } # end of loop
})

#Normal simulation methods

test_that("continuous simulation method matches a reference", {
 nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
 means.example <- c(22, 21, 21.5)
 sigma_sq.example <- c(1, 1, 0.9)
 sigma_b_sq.example <- c(0.1, 0.15, 0.1)
 multi.cps.normal.unbal <- cps.ma.normal(nsim = 100, nsubjects = nsubjects.example, 
                        means = means.example, sigma_sq = sigma_sq.example, 
                        sigma_b_sq = sigma_b_sq.example, alpha = 0.05,
                        quiet = FALSE, ICC=NULL, method = 'glmm', 
                        all.sim.data = FALSE,
                        seed = 123, cores = "all",
                        poor.fit.override = FALSE)
 expect_equal(round(multi.cps.normal.unbal[[1]][,1], 1), c(0.9, 0.3))
})

#doesn't pass, compare normal to t distributed rndom number generation
test_that("normal vs t-dist comparison", {
  q <- 10
  nc <- sample.int(200, q)
  ns <- sample.int(200, q)
  icc. <- runif(q, min=0.1, max=0.99)
  sig <- (runif(q, min=0.1, max=50))^2
  sigb <- createMissingVarianceParam(ICC = icc., sigma_b = sig)
  holder <- data.frame(nc, ns, icc., sig, sigb)
  same <- rep(NA, length=q)
  for (i in 1:q){
    multi.cps.norm <- cps.ma.normal(nsim = 200, narms = 2, 
                                      nclusters = nc[i], nsubjects = ns[i],
                                      means = c(0,1),
                                      tdist = FALSE,
                                      ICC = icc.[i],
                                      sigma_sq = sig[i], alpha = 0.05,
                                      quiet = FALSE, method = 'glmm',
                                      all.sim.data = FALSE,
                                      poor.fit.override = TRUE, cores="all")
    multi.cps.tdist <- cps.ma.normal(nsim = 200, narms = 2, 
                                      nclusters = nc[i], nsubjects = ns[i],
                                      means = c(0,1),
                                      tdist = TRUE,
                                      ICC = icc.[i],
                                      sigma_sq = sig[i], alpha = 0.05,
                                      quiet = FALSE, method = 'glmm',
                                      all.sim.data = FALSE,
                                      poor.fit.override = TRUE, cores="all")
    if(round(multi.cps.normal[[1]][,1], 1)==round(multi.cps.tdist[[1]][,1], 1)) {
      same[i] <- 1
    }  else {
      same[i] <- 0
    }
    #expect_equal(round(multi.cps.normal[[1]][,1], 1), round(as.numeric(analytic.mean)))
    print(paste("Interation", i, "of 10."))
  } # end of loop
})

 
#doesn't pass, analytic method is messed up
test_that("continuous simulation method matches the analytic method", {
  q <- 10
  nc <- sample.int(200, q)
  ns <- sample.int(200, q)
  icc. <- runif(q, min=0.1, max=0.99)
  sig <- (runif(q, min=0.1, max=50))^2
  sigb <- createMissingVarianceParam(ICC = icc., sigma_b = sig)
  holder <- data.frame(nc, ns, icc., sig, sigb)
  same <- rep(NA, length=q)
  for (i in 1:q){
    multi.cps.normal <- cps.ma.normal(nsim = 200, narms = 2, 
                                  nclusters = nc[i], nsubjects = ns[i],
                                  means = c(0,1),
                                  tdist = FALSE,
                                  ICC = icc.[i],
                                  sigma_sq = sig[i], alpha = 0.05,
                                  quiet = FALSE, method = 'glmm',
                                  all.sim.data = FALSE,
                                  poor.fit.override = TRUE, cores="all")
    twoarm.mean <- cps.normal(nsim = 200, nsubjects = ns[i], nclusters = nc[i], difference = 1,
                              ICC = icc.[i], sigma = sig[i], alpha = 0.05, method = 'glmm', 
                              quiet = FALSE, all.sim.data = FALSE)
    analytic.mean <- crtpwr.2mean(alpha = 0.05, power = NA, 
                              nclusters = nc[i], nsubjects = ns[i],
                              icc = icc.[i], 
                              vart = sig[i]+sigb[i], 
                              method = "weighted", 
                              tol = .Machine$double.eps^0.25, d=1)
  if(round(twoarm.mean$power[1], 1)==round(as.numeric(analytic.mean))) {
    same[i] <- 1
  }  else {
    same[i] <- 0
  }
  #expect_equal(round(multi.cps.normal[[1]][,1], 1), round(as.numeric(analytic.mean)))
  print(paste("Interation", i, "of 10."))
  } # end of loop
})

test_that("continuous simulation method matches the 2-arm simulation method", {
  nc <- sample.int(200, 10)
  ns <- sample.int(200, 10)
  icc. <- runif(10, min=0.1, max=0.99)
  sig <- runif(10, min=0.1, max=50)
  for (i in 1:10){
    multi.cps.normal <- cps.ma.normal(nsim = 100, narms = 2, 
                                      nclusters = nc[i], nsubjects = ns[i],
                                      means = c(0,1),
                                      ICC = icc.[i],
                                      sigma_sq = sig[i], alpha = 0.05,
                                      quiet = FALSE, method = 'glmm',
                                      all.sim.data = FALSE,
                                      poor.fit.override = TRUE, cores="all")
    twoarm.mean <- cps.normal(nsim = 100, nsubjects = ns[i], nclusters = nc[i], difference = 1,
                               ICC = icc.[i], sigma = sig[i], alpha = 0.05, method = 'glmm', 
                               quiet = FALSE, all.sim.data = FALSE)
    print(nc[i])
    print(ns[i])
    print(icc.[i])
    print(sig[i])
    expect_equal(round(multi.cps.normal[[1]][,1], 1), round(twoarm.mean$power[1], 1))
    print(paste("Interation", i, "of 10."))
  } # end of loop
})





## FIXME this is where I stopped
#compare balanced and unbalanced designs
 multi.cps.normal <- cps.ma.normal(nsim = 100, narms = 3, 
                                    nclusters = c (10,11,10), nsubjects = 100,
                                    means = c(21, 21, 21.4),
                                    sigma_sq = c(1,1,.9), 
                                    sigma_b_sq = c(.1,.15,.1), alpha = 0.05,
                                    quiet = FALSE, ICC=NULL, method = 'glmm',
                                    all.sim.data = FALSE, seed = 123,
                                    poor.fit.override = TRUE, cores="all")
 
 multi.cps.normal.simple <- cps.ma.normal(nsim = 100, narms = 3,
                                   nclusters = 10, nsubjects = 25, 
                                   means = c(22.1, 21, 22.5),
                                   sigma_sq = 1, 
                                   sigma_b_sq = 1, alpha = 0.05,
                                   quiet = FALSE, ICC=NULL, method = 'glmm',
                                   all.sim.data = FALSE, seed = 123,
                                   poor.fit.override = TRUE, cores="all")
 
 
 test_that("count simulation method matches the 2-arm simulation method", {
   nc <- sample.int(200, 10)
   ns <- sample.int(200, 10)
   sig <- runif(10, min=0.1, max=100)
   for (i in 1:1){
     count.ma <- cps.ma.count(nsim = 100, nsubjects = ns[i], narms = 2,
                                  nclusters = nc[i],
                                  counts = c(30, 70),
                                  sigma_b_sq = sig[i], alpha = 0.05,
                                  quiet = TRUE, method = 'glmm', 
                                  all.sim.data = FALSE, 
                                  multi.p.method="none",
                                  poor.fit.override = TRUE,
                                  cores="all")  
     count.sim = cps.count(nsim = 100, nsubjects = ns[i], nclusters = nc[i], c1 = 30,
                            c2 = 70, sigma_b = sig[i], family = 'poisson', analysis = 'poisson',
                            method = 'glmm', alpha = 0.05, quiet = FALSE, all.sim.data = TRUE)
     print(nc[i])
     print(ns[i])
     print(sig[i])
     expect_equal(as.numeric(round(count.ma[[1]][,1], 1)), as.numeric(round(count.sim$power[1], 1)))
     print(paste("Interation", i, "of 10."))
   } # end of loop
 })
