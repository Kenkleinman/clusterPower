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
#' @param nsubjects Number of subjects per cluster; accepts a single integer or a vector of 2 integers (if nsubjects differs between arms) (required). 
#' @param nclusters Number of clusters per treatment group; accepts a single integer or a vector of 2 integers (if nsubjects differs between arms) (required).
#' At least 2 of the following 3 arguments must be specified:
#' @param c1 Expected outcome count in non-treatment group
#' @param c2 Expected outcome count in treatment group
#' @param c.diff Expected difference in outcome count between groups, defined as c.diff = c1 - c2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal between groups. Accepts numeric
#' If between cluster variances differ between treatment groups, the following must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param family Distribution from which responses are simulated. Accepts Poisson ('poisson') or negative binomial ('neg.binom') (required); default = 'poisson'
#' @param analysis Family used for regression; currently only applicable for GLMM. Accepts c('poisson', 'neg.binom') (required); default = 'poisson'
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'
#' @param alpha Significance level. Default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default = FALSE.
#' @param all.sim.data Option to output list of all simulated datasets. Default = FALSE
#' 
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, distribution of 
#'   simulated data, and regression family
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Data frame containing families for distribution and analysis of simulated data
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting between-cluster variances for Treatment/Non-Treatment groups
#'   \item Vector containing expected counts and risk ratios based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value",
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)
#' }
#' @author Alexander R. Bogdan 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' 
#' @examples 
#' \dontrun{
#' count.sim = cps.count(nsim = 100, nsubjects = 20, nclusters = 36, c1 = 100,
#'                       c2 = 200, sigma_b = 1, family = 'poisson', analysis = 'poisson',
#'                       method = 'glmm', alpha = 0.05, quiet = FALSE, all.sim.data = TRUE)
#' }
#'
#' @export

# Define function
cps.count = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, c1 = NULL, c2 = NULL, 
                     c.diff = NULL, sigma_b = NULL, sigma_b2 = NULL, family = 'poisson', 
                     analysis = 'poisson', method = 'glmm', alpha = 0.05, quiet = FALSE, 
                     all.sim.data = FALSE, irgtt = FALSE, seed = NA){
  
  if (!is.na(seed)){
    set.seed(seed = seed)
  }
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  simulated.datasets = list()
  start.time = Sys.time()
  
  # Create progress bar
  prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta", 
                                         total = nsim, clear = FALSE, width = 100)
  prog.bar$tick(0)
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS
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
  # Set cluster sizes for treatment arm (if not already specified)
  if(length(nclusters) == 1){
    nclusters[2] = nclusters[1]
  }
  # Set sample sizes for each cluster (if not already specified)
  if(length(nsubjects) == 1){
    nsubjects = c(nsubjects, nsubjects)
  } 
  if(nclusters[1] == nclusters[2] && length(nsubjects) == nclusters[1]){
    nsubjects = rep(nsubjects, 2)
  }
  if(length(nclusters) == 2 && length(nsubjects) != 1 && length(nsubjects) != length(nclusters)){
    stop("A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS")
  }
  
  # Validate SIGMA_B, SIGMA_B2, ALPHA
  if (irgtt==FALSE){
    min0.warning = " must be a numeric value greater than 0"
    if(!is.numeric(sigma_b) || sigma_b <= 0){
      stop("SIGMA_B", min0.warning)
    }
    if(!is.null(sigma_b2) && sigma_b2 <= 0){
      stop("SIGMA_B2", min0.warning)
    }
  }
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
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
  
  # Validate FAMILY, ANALYSIS, METHOD, QUIET
  if(!is.element(family, c('poisson', 'neg.binom'))){
    stop("FAMILY must be either 'poisson' (Poisson distribution) 
         or 'neg.binom'(Negative binomial distribution)")
  }
  if(!is.element(analysis, c('poisson', 'neg.binom'))){
    stop("ANALYSIS must be either 'poisson' (Poisson regression) 
         or 'neg.binom'(Negative binomial regression)")
  }
  if(!is.element(method, c('glmm', 'gee'))){
    stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model) 
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  
  # Calculate inputs & variance parameters
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
  
  # Create indicators for treatment group & cluster
  trt = c(rep(0, length.out = nsubjects[1]*nclusters[1]), 
          #rep(1, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
          rep(1, length.out = nsubjects[2]*nclusters[2]))
  clust = c(rep(1:nclusters[1], each = nsubjects[1]), rep((nclusters[1] + 1):(nclusters[1] + nclusters[2]), 
                                                          each = nsubjects[2]))
  #clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))

  # Create simulation loop
  for(i in 1:nsim){
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b[2]))
    
    # Create non-treatment y-value
    #y0.intercept = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
    y0.intercept <- rep(randint.0, each = nsubjects[1])
    y0.linpred = y0.intercept + log(c1) 
    y0.prob = exp(y0.linpred)
    if(family == 'poisson'){
      y0 = stats::rpois(length(y0.prob), y0.prob)
    }
    if(family == 'neg.binom'){
      y0 = stats::rnbinom(length(y0.prob), size = 1, mu = y0.prob)
    }
      
    # Create treatment y-value
    #y1.intercept = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.intercept <- rep(randint.1, each = nsubjects[2])
    y1.linpred = y1.intercept + log(c2) #+ log((c1 / (1 - c1)) / (c2 / (1 - c2)))
    y1.prob = exp(y1.linpred)
    if(family == 'poisson'){
      y1 = stats::rpois(length(y1.prob), y1.prob)
    }
    if(family == 'neg.binom'){
      y1 = stats::rnbinom(length(y1.prob), size = 1, mu = y1.prob)
    }
    
    # Create single response vector
    y = c(y0, y1)
    
    # Create and store data for simulated dataset
    sim.dat = data.frame(y = as.integer(y), trt = as.factor(trt), clust = as.factor(clust))
    if(all.sim.data == TRUE){
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      if(irgtt == FALSE){
        if(analysis == 'poisson'){
          require("optimx")
          my.mod = lme4::glmer(y ~ trt + (1|clust), data = sim.dat, 
                               family = stats::poisson(link = 'log'))
        }
        if(analysis == 'neg.binom'){
          my.mod = lme4::glmer.nb(y ~ trt + (1|clust), data = sim.dat)
        }
        glmm.values = summary(my.mod)$coefficient
        est.vector = append(est.vector, glmm.values['trt1', 'Estimate'])
        se.vector = append(se.vector, glmm.values['trt1', 'Std. Error'])
        stat.vector = append(stat.vector, glmm.values['trt1', 'z value'])
        pval.vector = append(pval.vector, glmm.values['trt1', 'Pr(>|z|)'])  
      } else {
        if(analysis == 'poisson'){
         # require("optimx")
          my.mod <- lme4::glmer(y ~ trt + (0 + trt|clust), data = sim.dat, 
                               family = stats::poisson(link = 'log')#,                                
                              # control = lme4::glmerControl(optimizer = "optimx", 
          #                                                  optCtrl = list(method= "nlminb"))
          )
        }
        if(analysis == 'neg.binom'){
          my.mod = lme4::glmer.nb(y ~ trt + (0 + trt|clust), data = sim.dat)
        }        
      }
      glmm.values = summary(my.mod)$coefficient
      est.vector = append(est.vector, glmm.values['trt1', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt1', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt1', 'z value'])
      pval.vector = append(pval.vector, glmm.values['trt1', 'Pr(>|z|)'])  
    }
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y ~ trt, data = sim.dat,
                               family = stats::poisson(link = 'log'), 
                               id = clust, corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
    }
    
    # Update simulation progress information
    if(quiet == FALSE){
      if(i == 1){
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        message(paste0('Begin simulations :: Start Time: ', Sys.time(), ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
      }
      # Iterate progress bar
      prog.bar$update(i / nsim)
      Sys.sleep(1/100)
      
      if(i == nsim){
        message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
      } 
    }
  }
  
  ## Output objects
  # Create object containing summary statement
  if (irgtt==FALSE){
    summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, 
                             " Simulations: Simple Design, Count Outcome\nData Simulated from ", 
                             switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'), 
                             " distribution\nAnalyzed using ", 
                             switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'), 
                             " regression")
  } else {
    summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, 
                             " Simulations: IRGTT Design, Count Outcome\nData Simulated from ", 
                             switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'), 
                             " distribution\nAnalyzed using ", 
                             switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'), 
                             " regression")    
  }
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model', 
                       gee = 'Generalized Estimating Equation')
  
  # Store simulation output in data frame
  cps.model.est = data.frame(Estimate = as.vector(unlist(est.vector)),
                           Std.err = as.vector(unlist(se.vector)),
                           Test.statistic = as.vector(unlist(stat.vector)),
                           p.value = as.vector(unlist(pval.vector)))
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  pval.power = sum(cps.model.est[, 'sig.val']) / nrow(cps.model.est)
  power.parms = data.frame(power = round(pval.power, 3),
                           lower.95.ci = round(pval.power - abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           upper.95.ci = round(pval.power + abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  
  # Create object containing inputs
  c1.c2.rr = round(exp(log(c1) - log(c2)), 3)
  c2.c1.rr = round(exp(log(c2) - log(c1)), 3)
  inputs = t(data.frame('Non.Treatment' = c("count" = c1, "risk.ratio" = c1.c2.rr), 
                        'Treatment' = c("count" = c2, 'risk.ratio' = c2.c1.rr), 
                        'Difference' = c("count" = c.diff, 'risk.ratio' = c2.c1.rr - c1.c2.rr)))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Non.Treatment' = rep (nsubjects[1], length.out = nclusters[1]), 
                       'Treatment' = rep(nsubjects[2], length.out = nclusters[2]))
  
  # Create object containing number of clusters
  n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), 
                            "Treatment" = c("n.clust" = nclusters[2])))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame('Non.Treatment' = c('sigma_b' = sigma_b[1]), 
                           'Treatment' = c('sigma_b' = sigma_b[2])))
  
  # Create object containing FAMILY & REGRESSION parameters
  dist.parms = rbind('Family:' = paste0(switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'), ' distribution'), 
                     'Analysis:' = paste0(switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'), ' distribution'))
  colnames(dist.parms) = "Data Simuation & Analysis Parameters"
  
  
  # Create list containing all output (class 'crtpwr') and return
  complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = long.method, 
                                   "dist.parms" = dist.parms, "alpha" = alpha, "cluster.sizes" = cluster.sizes, 
                                   "n.clusters" = n.clusters, "variance.parms" = var.parms, "inputs" = inputs, 
                                   "model.estimates" = cps.model.est, "sim.data" = simulated.datasets), 
                              class = 'crtpwr')
  
  return(complete.output)
  }

