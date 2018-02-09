#' Power simulations for cluster-randomized trials: Simple Designs.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatments, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance; significance level, analytic method, and whether 
#' or not progress updates should be displayed while the function is running.
#' 
#'   
#' The \code{period.effect} parameter needs to be specified on the "link function scale". Meaning that if the
#' average baseline risk for a Poisson model is 4/1000, then the \code{period.effect} should be specified as
#' log(.004). Similarly, the baseline risk in a logistic model should be specified on the logit scale. The
#' period effect can have length of 1, in which case it is treated as the average period effect across all
#' periods, or it can have length equal to n.periods, in which case it is assumed that the investigator is
#' specifying exact period effects s/he wishes to simulate.
#' 
#' For the Poisson simulations, at risk time is computed for each individual in the simulation. If
#' \code{at.risk.time} is specified as a numeric vector of length 1, then the given number is the constant atrisk
#' time which every individual is assumed to have. If \code{length(at.risk.time)==2}, the values are taken
#' as the mean and size parameters of a negative binomial distribution (used as \code{mu} and \code{size} in the
#' \code{rnbinom()} function) from which an at-risk time is drawn for each individual. Specifically, the at risk
#' times are drawn as \code{at.risk.time = 1 + rnbinom(1, size=at.risk.params[2], mu=at.risk.params[1])}.
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param n Number of subjects per cluster; accepts integer (required). Currently in development to
#' accept vectors of treatment-specific cluster numbers and distributional forms. 
#' @param m Number of clusters per treatment group; accepts integer (required). Currently in development to accept 
#' a vector of length 2 to specify unique numbers of clusters by treatment group
#' @param difference Expected absolute treatment effect; accepts numeric (required).
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE); accepts c('glmm', 'gee') (required).
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default is FALSE.
#' At least 2 of the following must be specified:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma_w Within-cluster variance; accepts numeric
#' @param sigma_b Between-cluster variance; accepts numeric
#' If clusters differ between treatment groups, at least 2 of the following 
#' must be specified:
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma_w2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' 
#' @return A list with the following components
#' \describe{
#'   \item{sim.data}{Data frame with columns "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (t-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", "is.signif" (Is p-value less than alpha?)}
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
#' }
#' 
#' @author Alexander R. Bogdan
#' 
#' @examples 
#' \dontrun{
#' my.normal.sim = cps.normal(nsim = 100, n = 50, m = 6, difference = 30, ICC = 0.2, sigma_w = 100,
#'                     alpha = 0.05, method = 'glmm', quiet = FALSE)
#' }
#'
#' @export


# Define function
cps.normal = function(nsim = NULL, n = NULL, m = NULL, difference = NULL,
                      ICC = NULL, sigma_w = NULL, sigma_b = NULL, 
                      ICC2 = NULL, sigma_w2 = NULL, sigma_b2 = NULL, 
                      alpha = 0.05, method, quiet = FALSE){
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  start.time = Sys.time()
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, N, M, DIFFERENCE, ALPHA
  sim.data.arg.list = list(nsim, n, m, difference)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if(sum(sim.data.args) > 0){
    stop("NSIM, N, M & DIFFERENCE must all be specified. Please review your input values.")
  }
  min1.warning = " must be an integer greater than or equal to 1"
  if(!is.wholenumber(nsim) || nsim < 1){
    stop(paste0("NSIM", min1.warning))
  }
  if(!is.wholenumber(n) || n < 1){
    stop(paste0("N", min1.warning))
  }
  if(!is.wholenumber(m) || m < 1){
    stop(paste0("M", min1.warning))
  }
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(difference) || difference < 0){
    stop("DIFFERENCE", min0.warning)
  }
  if(!is.numeric(alpha) || alpha < 0){
    stop("ALPHA", min0.warning)
  } else if(alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate ICC, sigma_w, sigma_b, ICC2, sigma_w2, sigma_b2
  parm1.arg.list = list(ICC, sigma_w, sigma_b)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma_w, sigma_b")
  }
  parm2.arg.list = list(ICC2, sigma_w2, sigma_b2)
  parm2.args = unlist(lapply(parm2.arg.list, is.null))
  if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
    stop("At least two of the following terms must be provided to simulate treatment-specific
         variances: ICC2, sigma_w2, sigma_b2")
  }
  
  # Validate METHOD, QUIET
  if(!is.element(method, c('glmm', 'gee'))){
    stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model) 
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  
  # Create simulation loop
  for(i in 1:nsim){
    ## Create simulated response variable
    # Single set of parameters
    if(!is.null(ICC)){
      icc_jitter = ICC * sqrt((n - m[1]) / stats::rchisq(1, df = n - m[1]))
    }
    if(!is.null(c(ICC, sigma_w)) && is.null(sigma_b)){
      sigma_b = icc_jitter * sigma_w / (1 - icc_jitter)
    }
    if(!is.null(c(ICC, sigma_b)) && is.null(sigma_w)){
      sigma_w = sigma_b / icc_jitter - sigma_b
    }

    # Second set of parameters
    if(!is.null(ICC2)){
      icc2_jitter = ICC2 * sqrt((n - m[2]) / stats::rchisq(1, df = n - m[2]))
    }
    if(!is.null(c(ICC2, sigma_w2)) && is.null(sigma_b2)){
      sigma_b2 = icc2_jitter * sigma_w2 / (1 - icc2_jitter)
    }
    if(!is.null(c(ICC2, sigma_b2)) && is.null(sigma_w2)){
      sigma_w2 = sigma_b2 / icc2_jitter - sigma_b2
    }
    
    # Generate simulated data with treatment-specific cluster variances
    if(!is.null(c(sigma_w2, sigma_b2))){
      # Generate simulated response values with sigma_w & sigma_b
      y.vals = NULL
      sim.base = stats::rchisq(m, df = m - 1)
      clust.means = append(y.vals, apply(as.matrix(sim.base), 1, function(x) stats::rchisq(n, df = n - 1, ncp = x)))
      within.clust = stats::rt(n*m, df = n - 1, ncp = sqrt(sigma_w))
      between.clust = stats::rt(m, df = m - 1, ncp = sqrt(sigma_b))
      y = clust.means + within.clust + between.clust
      
      # Generate simulated response values with sigma_w2 & sigma_b2
      y.vals = NULL
      sim.base = stats::rchisq(m, df = m - 1)
      clust.means = append(y.vals, apply(as.matrix(sim.base), 1, function(x) stats::rchisq(n, df = n - 1, ncp = x)))
      within.clust = stats::rt(n*m, df = n - 1, ncp = sqrt(sigma_w2))
      between.clust = stats::rt(m, df = m - 1, ncp = sqrt(sigma_b2))
      y.trt = clust.means + within.clust + between.clust + difference
      y = append(y, y.trt)
    }
    
    # Generate simulated data with global cluster variances
    if(is.null(c(sigma_w2, sigma_b2))){
      # Generate simulated response values with sigma_w & sigma_b
      y.vals = NULL
      sim.base = stats::rchisq(m*2, df = 2*m - 2)
      clust.means = append(y.vals, apply(as.matrix(sim.base), 1, function(x) stats::rchisq(n, df = n - 1, ncp = x)))
      within.clust = stats::rt(n*m*2, df = n - 1, ncp = sqrt(sigma_w))
      between.clust = stats::rt(m*2, df = 2*m - 2, ncp = sqrt(sigma_b))
      y = clust.means + within.clust + between.clust
      y[(1+n*m):(n*m*2)] = y[(1+n*m):(n*m*2)] + difference
    }
  
    # Create treatment and cluster indicator variables
    trt = c(rep(0, n*m), rep(1, n*m))
    clust = c(rep(1:m, n), rep((m+1):(m*2), n))
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(y.resp = y, y = y, trt = trt, clust = clust)
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y.resp ~ trt + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = m * 2 - 2)
      est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt', 't value'])
      pval.vector = append(pval.vector, p.val)
    }
    # Fit GEE (geeglm)
    if(method == 'gee'){
      my.mod = geepack::geeglm(y.resp ~ trt, data = sim.dat, 
                               id = clust, corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
    }
    
    if(quiet == FALSE){
      if(i == 1){
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        message(paste0('Begin simulations :: Start Time: ', Sys.time(), ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
      }
      else if(i == nsim){
        message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
      } 
      else if(i %% 10 == 0){
        time.est = avg.iter.time * (nsim - i) / 60 
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        min.est = ifelse(min.est == 0, '<1', min.est)
        message(paste0('Progress: ', i / nsim * 100, '% complete :: Estimated time remaining: ', hr.est, 'Hr:', min.est, 'Min'))
      }
    }
  }
  cps.sim.dat = data.frame(estimates = as.vector(unlist(est.vector)), 
                           stderrs = as.vector(unlist(se.vector)),
                           test.stat = as.vector(unlist(stat.vector)),
                           pvals = as.vector(unlist(pval.vector)))
  cps.sim.dat[, 'sig.vals'] = ifelse(cps.sim.dat[, 'pvals'] < alpha, 1, 0)
  pval.power = sum(cps.sim.dat[, 'sig.vals']) / length(cps.sim.dat[, 'sig.vals'])
  power.parms = data.frame(power = round(pval.power, 3), 
                           lower.95.ci = round(pval.power - abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3), 
                           upper.95.ci = round(pval.power + abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  complete.output = list("sim.data" = cps.sim.dat, "power" = power.parms)
  return(complete.output)
  }

# Test function
my.sim = cps.normal(nsim = 100, n = 500, m = 30, difference = 3, 
                    ICC = 0.2, sigma_w = 100, ICC2 = 0.2, sigma_w2 = 100, 
                    alpha = 0.05, method = 'gee', quiet = FALSE)
my.sim$power

