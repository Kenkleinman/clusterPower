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
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters per treatment group; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified:
#' @param c1 Expected outcome count in non-treatment group
#' @param c2 Expected outcome count in treatment group
#' @param c.diff Expected difference in outcome count between groups, defined as c.diff = c1 - c2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal between groups. Accepts numeric
#' If between cluster variances differ between treatment groups, the following must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @alpha Significance level for power estimation, accepts value between 0 - 1; default = 0.05
#' @param family Distribution from which responses are simulated. Accepts c('poisson', 'neg.binom') (required); default = 'poisson'
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'
#' @param alpha Significance level. Default = 0.05.
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
#' my.count.sim = cps.count(nsim = 100, nsubjects = 50, nclusters = 6, c1 = 100/1000, c2 = 25/1000, sigma_b = 100,
#'                     family = 'poisson', method = 'glmm', alpha = 0.05, quiet = FALSE)
#' }
#'
#' @export

# Define function
cps.count = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, c1 = NULL, c2 = NULL, 
                     c.diff = NULL, sigma_b = NULL, sigma_b2 = NULL, 
                     family = 'poisson', method = 'glmm', alpha = 0.05, quiet = FALSE){
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  start.time = Sys.time()
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, N, M, SIGMA_B, ALPHA
  sim.data.arg.list = list(nsim, nclusters, nsubjects, sigma_b)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if(sum(sim.data.args) > 0){
    stop("NSIM, NSUBJECTS, NCLUSTERS & SIGMA_B must all be specified. Please review your input values.")
  }
  min1.warning = " must be an integer greater than or equal to 1"
  if(!is.wholenumber(nsim) || nsim < 1){
    stop(paste0("NSIM", min1.warning))
  }
  if(!is.wholenumber(nclusters) || nclusters < 1){
    stop(paste0("NCLUSTERS", min1.warning))
  }
  if(!is.wholenumber(nsubjects) || nsubjects < 1){
    stop(paste0("NSUBJECTS", min1.warning))
  }
  if(length(nclusters) > 2){
    stop("NCLUSTERS can only be a vector of length 1 (equal # of clusters per group) or 2 (unequal # of clusters per group)")
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
  if(length(nclusters) == 1){
    nclusters[2] = nclusters[1]
  }
  
  # Set sample sizes for each cluster (if not already specified)
  if(length(nsubjects) == 1){
    nsubjects[1:sum(nclusters)] = nsubjects
  } 
  if(nclusters[1] == nclusters[2] && length(nsubjects) == nclusters[1]){
    nsubjects = rep(nsubjects, 2)
  }
  if(length(nclusters) == 2 && length(nsubjects) != 1 && length(nsubjects) != sum(nclusters)){
    stop("A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS")
  }
  
  # Validate C1, C2, C.DIFF
  parm1.arg.list = list(c1, c2, c.diff)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: C1, C2, C.DIFF")
  }
  if(sum(parm1.args) == 0 && c.diff != abs(c1 - c2)){
    stop("At least one of the following terms has be misspecified: C1, C2, C.DIFF")
  }
  
  # Validate FAMILY, METHOD, QUIET
  if(!is.element(family, c('poisson', 'neg.binom'))){
    stop("FAMILY must be either 'poisson' (Poisson distribution) 
         or 'neg.binom'(Negative binomial distribution)")
  }
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
  if(is.null(c.diff)){
    c.diff = c1 - c2
  }
  if(is.null(sigma_b2)){
    sigma_b[2] = sigma_b
  }else{
    sigma_b[2] = sigma_b2
  }
  
  # Create simulation loop
  for(i in 1:nsim){
    # Create indicators for treatment group & cluster
    trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])), 
            rep(1, length.out = sum(nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])])))
    clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))
    
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b[2]))
    
    # Create non-treatment y-value
    y0.intercept = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
    y0.linpred = y0.intercept + log(c1)
    y0.prob = exp(y0.linpred)
    if(family == 'poisson'){
      y0 = rpois(length(y0.prob), y0.prob)
    }
    if(family == 'neg.binom'){
      y0 = rnbinom(length(y0.prob), size = 1, mu = y0.prob)
    }
      
    # Create treatment y-value
    y1.intercept = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1]+x])))
    y1.linpred = y1.intercept + log(c2) + log((c1 / (1 - c1)) / (c2 / (1 - c2)))
    y1.prob = exp(y1.linpred)
    if(family == 'poisson'){
      y1 = rpois(length(y1.prob), y1.prob)
    }
    if(family == 'neg.binom'){
      y1 = rnbinom(length(y1.prob), size = 1, mu = y1.prob)
    }
    
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
  # Create object containing summary statement
  summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, 
                           " Simulations: Count Outcome\nData Simulated from ", 
                           switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'), 
                           " distribution")
  
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
  
  # Create object containing expected difference
  difference = paste0("Expected Count Difference = ", abs(c.diff))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Group 1 (Non-Treatment)' = nsubjects[1:nclusters[1]], 
                       'Group 2 (Treatment)' = nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), "Treatment" = c("n.clust" = nclusters[2])))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame('Group.1.Non.Treatment' = c('sigma_b' = sigma_b[1]), 
                           'Group.2.Treatment' = c('sigma_b' = sigma_b[2])))
  
  # Create list containing all output and return
  complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = method, "alpha" = alpha,
                                   "cluster.sizes" = cluster.sizes, "n.clusters" = n.clusters, "variance.parms" = var.parms, 
                                   "difference" = difference, "sim.data" = cps.sim.dat), class = 'crtpwr')
  
  return(complete.output)
  }
