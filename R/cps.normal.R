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
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param m Number of subjects per cluster; accepts integer (required). Currently in development to
#' accept vectors of treatment-specific cluster numbers and distributional forms. 
#' @param n Number of clusters per treatment group; accepts single integer (required) or 
#' vector of length 2 for unequal number of clusters per treatment group 
#' @param difference Expected absolute treatment effect; accepts numeric (required).
#' At least 2 of the following must be specified:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma_w Within-cluster variance; accepts numeric
#' @param sigma_b Between-cluster variance; accepts numeric
#' If clusters differ between treatment groups, at least 2 of the following 
#' must be specified:
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma_w2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level. Default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default is FALSE.

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
#' my.normal.sim = cps.normal(nsim = 100, m = 50, n = 30, difference = 30, ICC = 0.2, sigma_w = 100,
#'                     alpha = 0.05, method = 'glmm', quiet = FALSE)
#' }
#'
#' @export


cps.normal = function(nsim = NULL, m = NULL, n = NULL, difference = NULL,
                      ICC = NULL, sigma = NULL, sigma_b = NULL,
                      ICC2 = NULL, sigma2 = NULL, sigma_b2 = NULL,
                      alpha = 0.05, method = 'glmm', quiet = FALSE){
  
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  
  # Set start.time for progress iterator
  start.time = Sys.time()
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, N, M
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
  if(length(n) > 2){
    stop("N can only be a vector of length 1 (equal # of clusters per group) or 2 (unequal # of clusters per group)")
  }
  
  # Set cluster sizes for treatment arm (if not already specified)
  if(length(n) == 1){
    n[2] = n[1]
  }
  
  # Set sample sizes for each cluster (if not already specified)
  if(length(m) == 1){
    m[1:sum(n)] = m
  } 
  if(n[1] == n[2] && length(m) == n[1]){
    m = rep(m, 2)
  }
  if(length(n) == 2 && length(m) != 1 && length(m) != sum(n)){
    stop("A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for M")
  }
  
  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(difference) || difference < 0){
    stop("DIFFERENCE", min0.warning)
  }
  
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate ICC, sigma, sigma_b, ICC2, sigma2, sigma_b2
  parm1.arg.list = list(ICC, sigma, sigma_b)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma, sigma_b")
  }
  parm2.arg.list = list(ICC2, sigma2, sigma_b2)
  parm2.args = unlist(lapply(parm2.arg.list, is.null))
  if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
    stop("At least two of the following terms must be provided to simulate treatment-specific
         variances: ICC2, sigma2, sigma_b2")
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
    # Single set of cluster parameters
    if(!is.null(c(ICC, sigma)) && is.null(sigma_b)){
      sigma_b = ICC * sigma / (1 - ICC)
    }
    if(!is.null(c(ICC, sigma_b)) && is.null(sigma)){
      sigma = sigma_b / ICC - sigma_b
    }
    if(!is.null(c(sigma, sigma_b)) && is.null(ICC)){
      ICC = sigma_b / (sigma_b + sigma)
    }
    
    # Second set of cluster parameters
    if(!is.null(c(ICC2, sigma2)) && is.null(sigma_b2)){
      sigma_b2 = ICC2 * sigma2 / (1 - ICC2)
    }
    if(!is.null(c(ICC2, sigma_b2)) && is.null(sigma2)){
      sigma2 = sigma_b2 / ICC2 - sigma_b2
    }
    if(!is.null(c(sigma2, sigma_b2)) && is.null(ICC2)){
      ICC2 = sigma_b2 / (sigma_b2 + sigma2)
    }
    
    # Set within/between cluster variance & ICC for treatment group (if not already specified)
    sigma[2] = ifelse(!is.null(sigma2), sigma2, sigma[1])
    sigma_b[2] = ifelse(!is.null(sigma_b2), sigma_b2, sigma_b[1])
    ICC[2] = ifelse(!is.null(ICC2), ICC2, ICC[1])
    
    # Generate simulated data
    # Create indicators for treatment group & cluster
    trt = c(rep(0, length.out = sum(m[1:n[1]])), rep(1, length.out = sum(m[(n[1]+1):(n[1]+n[2])])))
    clust = unlist(lapply(1:sum(n), function(x) rep(x, length.out = m[x])))
    
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(n[1], mean = 0, sd = sqrt(sigma_b[1]))
    randint.1 = stats::rnorm(n[2], mean = 0, sd = sqrt(sigma_b[2]))
    
    # Create non-treatment y-value
    y0.bclust = unlist(lapply(1:n[1], function(x) rep(randint.0[x], length.out = m[x])))
    y0.wclust = unlist(lapply(m[1:n[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[1]))))
    y.0 = y0.bclust + y0.wclust
    
    # Create treatment y-value
    y1.bclust = unlist(lapply(1:n[2], function(x) rep(randint.1[x], length.out = m[n[1]+x])))
    y1.wclust = unlist(lapply(m[(n[1]+1):(n[1]+n[2])], function(x) stats::rnorm(x, mean = difference, sd = sqrt(sigma[2]))))
    y.1 = y1.bclust + y1.wclust
    
    # Create single response vector
    y = c(y.0,y.1)
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(y.resp = y, trt = trt, clust = clust)
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y.resp ~ trt + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = sum(n) - 1)
      est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt', 't value'])
      pval.vector = append(pval.vector, p.val)
    }
    
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y.resp ~ trt, data = sim.dat,
                               id = clust, corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
    }
    
    # Create & update progress bar
    if(quiet == FALSE){
      if(i == 1){
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        message(paste0('Begin simulations :: Start Time: ', Sys.time(), ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
      }
      else if(i == nsim){
        total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        hr.est = total.est %/% 3600
        min.est = total.est %/% 60
        sec.est = round(total.est %% 60, 0)
        message(paste0("Simulations Complete! Time Completed: ", Sys.time(), "\nTotal Runtime: ", hr.est, 'Hr:', min.est, 'Min:', sec.est, 'Sec'))
      }
      else if(as.numeric(avg.iter.time) > 10 && i %% 5 == 0){
        time.est = avg.iter.time * (nsim - i) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        min.est = ifelse(min.est == 0, '<1', min.est)
        message(paste0('Progress: ', i / nsim * 100, '% complete :: Estimated time remaining: ', hr.est, 'Hr:', min.est, 'Min'))
      }
      else if(as.numeric(avg.iter.time) < 1 && i %% round(nsim/10, 0) == 0){
        time.est = avg.iter.time * (nsim - i) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        min.est = ifelse(min.est == 0, '<1', min.est)
        message(paste0('Progress: ', i / nsim * 100, '% complete :: Estimated time remaining: ', hr.est, 'Hr:', min.est, 'Min'))
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
  # Create object containing summary statement
  summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, " Simulations")
  
  # Store simulation output in data frame
  cps.sim.dat = data.frame(estimates = as.vector(unlist(est.vector)),
                           stderrs = as.vector(unlist(se.vector)),
                           test.stat = as.vector(unlist(stat.vector)),
                           pvals = as.vector(unlist(pval.vector)))
  cps.sim.dat[, 'sig.vals'] = ifelse(cps.sim.dat[, 'pvals'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  pval.power = sum(cps.sim.dat[, 'sig.vals']) / nrow(cps.sim.dat)
  power.parms = data.frame(power = round(pval.power, 3),
                           lower.95.ci = round(pval.power - abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           upper.95.ci = round(pval.power + abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Group 1 (Non-Treatment)' = m[1:n[1]], 'Group 2 (Treatment)' = m[(n[1]+1):(n[1]+n[2])])
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame('Group.1.Non.Treatment' = c('ICC' = ICC[1], 'sigma' = sigma[1], 'sigma_b' = sigma_b[1]), 
                           'Group.2.Treatment' = c('ICC' = ICC[2], 'sigma' = sigma[2], 'sigma_b' = sigma_b[2])))
    
  # Create list containing all output and return
  complete.output = list("overview" = summary.message, "sim.data" = cps.sim.dat, "power" = power.parms, "cluster.sizes" = cluster.sizes, 
                         "variance.parms" = var.parms, "alpha" = alpha)
  return(complete.output)
  }
