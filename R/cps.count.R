#' Power simulations for cluster-randomized trials: Simple Designs, Count Outcome.
#'
#' This function utilizes iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation for count outcomes.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, between-cluster variance, 
#' two of the following: expected count in non-treatment group, expected count 
#' in treatment group, difference in counts between groups; significance level, 
#' analytic method, and whether or not progress updates should be displayed 
#' while the function is running.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param m Number of subjects per cluster; accepts integer (required). 
#' @param n Number of clusters per treatment group; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified:
#' @param c1 Expected outcome count in non-treatment group
#' @param c2 Expected outcome count in treatment group
#' @param c.diff Expected difference in outcome count between groups, defined as c.diff = c1 - c2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal between groups. Accepts numeric
#' If between cluster variances differ between treatment groups, the following must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @alpha Significance level for power estimation, accepts value between 0 - 1; default = 0.05
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE); accepts c('glmm', 'gee') (required).
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default = FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{sim.data}{Data frame with columns "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
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
#' my.count.sim = cps.count(nsim = 100, n = 50, m = 6, c1 = 100, c2 = 25, sigma_b = 100,
#'                     alpha = 0.05, method = 'glmm', quiet = FALSE)
#' }
#'
#' @export

# Define function
cps.count = function(nsim = NULL, m = NULL, n = NULL, c1 = NULL, c2 = NULL, 
                     c.diff = NULL, sigma_b = NULL, sigma_b2 = NULL, 
                     alpha = 0.05, method, quiet = FALSE){
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  start.time = Sys.time()
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, N, M, SIGMA_B, ALPHA
  sim.data.arg.list = list(nsim, n, m, sigma_b)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if(sum(sim.data.args) > 0){
    stop("NSIM, N, M & SIGMA_B must all be specified. Please review your input values.")
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
  
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(sigma_b) || sigma_b <= 0){
    stop("SIGMA_B", min0.warning)
  }
  if(!is.null(sigma_b2) && sigma_b2 <= 0){
    stop("SIGMA_B2", min0.warning)
  }
  
  if(!is.numeric(alpha) || alpha < 0){
    stop("ALPHA", min0.warning)
  } else if(alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
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
  
  # Validate C1, C2, C.DIFF
  parm1.arg.list = list(c1, c2, c.diff)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: C1, C2, C.DIFF")
  }
  
  # Validate METHOD, QUIET
  if(!is.element(method, c('glmm', 'gee'))){
    stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model) 
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  
  # Simulation parameters
  if(is.null(c1)){
    c1 = abs(c.diff - c2)
  }
  if(is.null(c2)){
    c2 = abs(c1 - c.diff)
  }
  #if(is.null(c.diff)){
  #  c.diff = c1 - c2
  #}
  if(is.null(sigma_b2)){
    sigma_b[2] = sigma_b
  }else{
    sigma_b[2] = sigma_b2
  }
  
  # Set beta
  beta = 0.5
  
  # Create simulation loop
  for(i in 1:nsim){
    # Generate simulated data
    # Create indicators for treatment group & cluster
    trt = c(rep(0, length.out = sum(m[1:n[1]])), rep(1, length.out = sum(m[(n[1]+1):(n[1]+n[2])])))
    clust = unlist(lapply(1:sum(n), function(x) rep(x, length.out = m[x])))
    
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(n[1], mean = 0, sd = sqrt(sigma_b[1]))
    randint.1 = stats::rnorm(n[2], mean = 0, sd = sqrt(sigma_b[2]))
    
    # Create non-treatment y-value
    y0.intercept = unlist(lapply(1:n[1], function(x) rep(randint.0[x], length.out = m[x])))
    y0.test = unlist(lapply(m[1:n[1]], function(x) stats::rnorm(x, c1)))
    y0.linpred = y0.intercept + y0.test * beta
    y0.prob = exp(y0.linpred)
    y0 = rpois(length(y0.prob), y0.prob)
    
    # Create treatment y-value
    y1.intercept = unlist(lapply(1:n[2], function(x) rep(randint.1[x], length.out = m[n[1]+x])))
    y1.test = unlist(lapply(m[(n[1]+1):(n[1]+n[2])], function(x) stats::rnorm(x, c2)))
    y1.linpred = y1.intercept + y1.test * beta
    y1.prob = exp(y1.linpred)
    y1 = rpois(length(y1.prob), y1.prob)
    
    # Create single response vector
    y = c(y0,y1)
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(y.resp = y, trt = trt, clust = clust)
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::glmer(y.resp ~ trt + (1|clust), data = sim.dat, family = poisson(link = 'log'))
      #suppressWarnings(lme4::glmer(y.resp ~ trt + (1|clust), data = sim.dat, family = poisson(link = 'log')))
      glmm.values = summary(my.mod)$coefficient
      est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt', 'z value'])
      pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|z|)'])
    }
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y.resp ~ trt, data = sim.dat,
                               family = poisson(link = 'log'), 
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
my.sim = cps.count(nsim = 100, m = 50, n = 5, c1 = 60/100, c2 = 25/100, sigma_b = 25, alpha=0.05, method='glmm')
my.sim$power

