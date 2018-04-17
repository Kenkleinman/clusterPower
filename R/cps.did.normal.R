#' Power simulations for cluster-randomized trials: Difference in Difference Design, Continuous Outcome.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for difference in difference cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation for difference in difference (DID) cluster-randomized controlled trial.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatment arms, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per treatment group; accepts either a scalar (equal cluster sizes, both groups), 
#' a vector of length two (equal cluster sizes within groups), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within groups) (required).
#' accept vectors of treatment-specific cluster numbers and distributional forms.
#' @param nclusters Number of clusters per group; accepts single integer or vector of length 2 for unequal number 
#' of clusters per treatment group (required)
#' @param difference Expected absolute difference in treatment effect between time points; accepts numeric (required).
#' The following arguments can accept either a scalar (indicating constant variance between time points) 
#' or a vector of length 2 (indicating different variances at each time point). 
#' At least two of the following must be specified for the non-treatment group:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma_w Within-cluster variance; accepts numeric
#' @param sigma_b Between-cluster variance; accepts numeric
#' Similarly, at least two of the following must be specified for the treatment group:
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma_w2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level. Default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{overview}{Character string indicating total number of simulations and simulation type}
#'   \item{nsim}{Number of simulations}
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
#'   \item{method}{Analytic method used for power estimation}
#'   \item{alpha}{Significance level}
#'   \item{cluster.sizes}{Vector containing user-defined cluster sizes}
#'   \item{n.clusters}{Vector containing user-defined number of clusters}
#'   \item{variance.parms}{Data frame reporting ICC for Treatment/Non-Treatment groups at each time point}
#'   \item{inputs}{Vector containing expected difference between groups based on user inputs}
#'   \item{model.estimates}{Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)}
#'   \item{sim.data}{List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster), 
#'                   "period" (Indicator for time)}
#' }
#' 
#' @author Alexander R. Bogdan
#' 
#' @examples 
#' \dontrun{
#' my.did.rct = cps.did.normal(nsim = 100, nsubjects = 50, nclusters = 30, difference = 3, ICC = 0.2, sigma = 100,
#'                     alpha = 0.05, method = 'glmm', quiet = FALSE, all.sim.data = FALSE)
#' }
#'
#' @export


cps.did.normal = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, difference = NULL,
                      ICC = NULL, sigma = NULL, sigma_b = NULL,
                      ICC2 = NULL, sigma2 = NULL, sigma_b2 = NULL,
                      alpha = 0.05, method = 'glmm', quiet = FALSE,
                      all.sim.data = FALSE){
  
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  simulated.datasets = list()
  
  # Set start.time for progress iterator & initialize progress bar
  start.time = Sys.time()
  prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta", 
                                         total = nsim, clear = FALSE, width = 100)
  prog.bar$tick(0)
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Validate NSIM, N, M
  sim.data.arg.list = list(nsim, nclusters, nsubjects, difference)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if(sum(sim.data.args) > 0){
    stop("NSIM, NCLUSTERS, NSUBJECTS & DIFFERENCE must all be specified. Please review your input values.")
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
    nsubjects[1:sum(nclusters)] = nsubjects
  }
  if(length(nsubjects) == 2){
    nsubjects = c(rep(nsubjects[1], nclusters[1]), rep(nsubjects[2], nclusters[2]))
  }
  if(nclusters[1] == nclusters[2] && length(nsubjects) == nclusters[1]){
    nsubjects = rep(nsubjects, 2)
  }
  if(length(nclusters) == 2 && length(nsubjects) != 1 && length(nsubjects) != sum(nclusters)){
    stop("A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS")
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
  if(sum(parm1.args) == 0 && ICC != sigma_b / (sigma_b + sigma)){
    stop("At least one of the following terms has be misspecified: ICC, sigma, sigma_b")
  }
  parm2.arg.list = list(ICC2, sigma2, sigma_b2)
  parm2.args = unlist(lapply(parm2.arg.list, is.null))
  if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
    stop("At least two of the following terms must be provided to simulate treatment-specific
         variances: ICC2, sigma2, sigma_b2")
  }
  if(sum(parm2.args) == 0 && ICC2 != sigma_b2 / (sigma_b2 + sigma2)){
    stop("At least one of the following terms has be misspecified: ICC2, sigma2, sigma_b2")
  }
  
  # Validate METHOD, QUIET, ALL.SIM.DATA
  if(!is.element(method, c('glmm', 'gee'))){
    stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model)
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  if(!is.logical(all.sim.data)){
    stop("ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output")
  }
  
  ## Create variance parameters
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
  
  # Create indicators for time, treatment group & cluster
  period = rep(0:1, each = 2 * sum(nsubjects))
  trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])), 
          rep(1, length.out = sum(nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))
  
  #####################################################
  ##########     SIMULATION ENGINE          ###########
  #####################################################
  
  # Create simulation parameters
  # nsubjects = rep(50, 60)
  # nclusters = c(30, 30)
  # difference = 3
  # 
  # sigma0 = c(10, 10)
  # sigma1 = c(10, 10)
  # sigma_b0 = c(15, 15)
  # sigma_b1 = c(15, 15)
  # 
  # # Create indicators for time, treatment group & cluster
  # period = rep(0:1, each = sum(nsubjects))
  # trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])), 
  #         rep(1, length.out = sum(nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])])))
  # clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))
  # 
  #   ### Generate simulated data
  #   ## TIME == 0
  #   # Generate between-cluster effects for non-treatment and treatment
  #   randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b0[1]))
  #   randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b0[2]))
  #   
  #   # Create non-treatment y-value
  #   y0.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
  #   y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma0[1]))))
  #   y0.pre = y0.bclust + y0.wclust + rnorm(nsubjects[1:nclusters[1]])
  #   
  #   # Create treatment y-value
  #   y1.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
  #   y1.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])], 
  #                             function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma0[2]))))
  #   y1.pre = y1.bclust + y1.wclust + rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  #   
  #   ## TIME == 1
  #   # Generate between-cluster effects for non-treatment and treatment
  #   randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b1[1]))
  #   randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b1[2]))
  #   
  #   # Create non-treatment y-value
  #   y0.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
  #   y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma1[1]))))
  #   y0.post = y0.bclust + y0.wclust + rnorm(nsubjects[1:nclusters[1]])
  #   
  #   # Create treatment y-value
  #   y1.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
  #   y1.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])], 
  #                             function(x) stats::rnorm(x, mean = difference, sd = sqrt(sigma1[2]))))
  #   y1.post = y1.bclust + y1.wclust + rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  #   
  #   # Create single response vector
  #   y = c(y0.pre, y1.pre, y0.post, y1.post)
  #   
  #   # Create data frame for simulated dataset
  #   sim.dat = data.frame(y = y, trt = trt, clust = clust, period = period)
  #   
  #   library(dplyr)
  #   sim.dat %>%
  #     group_by(period, trt) %>%
  #     summarise(avg = mean(y))
  #   
  # 
  # #####################################################
  
  
  # Create simulation loop
  for(i in 1:nsim){
    ### Generate simulated data
    ## TIME == 0
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b[2]))
    
    # Create non-treatment y-value
    y0.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
    y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[1]))))
    y0.pre = y0.bclust + y0.wclust
    
    # Create treatment y-value
    y1.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1]+x])))
    y1.wclust = unlist(lapply(nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])], 
                              function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[2]))))
    y1.pre = y1.bclust + y1.wclust
    
    ## TIME == 1
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b2[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b2[2]))
    
    # Create non-treatment y-value
    y0.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
    y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[1]))))
    y0.post = y0.bclust + y0.wclust
    
    # Create treatment y-value
    y1.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1]+x])))
    y1.wclust = unlist(lapply(nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])], 
                              function(x) stats::rnorm(x, mean = difference, sd = sqrt(sigma[2]))))
    y1.post = y1.bclust + y1.wclust
    
    # Create single response vector
    y = c(y0.pre, y1.pre, y0.post, y1.post)
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(y = y, trt = trt, clust = clust, period = period)
    if(all.sim.data == TRUE){
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + period + trt:period + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = sum(nclusters) - 2)
      est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt', 't value'])
      pval.vector = append(pval.vector, p.val)
    }
    
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y ~ trt + period + trt:period, data = sim.dat,
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
        message(paste0('Begin simulations :: Start Time: ', Sys.time(), 
                       ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
      }
      
      # Iterate progress bar
      prog.bar$update(i / nsim)
      Sys.sleep(1/100)
      
      if(i == nsim){
        total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        hr.est = total.est %/% 3600
        min.est = total.est %/% 60
        sec.est = round(total.est %% 60, 0)
        message(paste0("Simulations Complete! Time Completed: ", Sys.time(), 
                       "\nTotal Runtime: ", hr.est, 'Hr:', min.est, 'Min:', sec.est, 'Sec'))
      }
    }
  }
  # Create object containing summary statement
  summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, " Simulations: Difference in Difference, Continuous Outcome")
  
  # Store simulation output in data frame
  cps.model.est = data.frame(Estimate = as.vector(unlist(est.vector)),
                             Std.err = as.vector(unlist(se.vector)),
                             Test.statistic = as.vector(unlist(stat.vector)),
                             p.value = as.vector(unlist(pval.vector)))
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  pval.power = sum(cps.model.est[, 'sig.val']) / nrow(cps.model.est)
  power.parms = data.frame(Power = round(pval.power, 3),
                           Lower.95.CI = round(pval.power - abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           Upper.95.CI = round(pval.power + abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Group 1 (Non-Treatment)' = nsubjects[1:nclusters[1]], 
                       'Group 2 (Treatment)' = nsubjects[(nclusters[1]+1):(nclusters[1]+nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), "Treatment" = c("n.clust" = nclusters[2])))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame('Non.Treatment' = c('ICC' = ICC[1], 'sigma' = sigma[1], 'sigma_b' = sigma_b[1]), 
                           'Treatment' = c('ICC' = ICC[2], 'sigma' = sigma[2], 'sigma_b' = sigma_b[2])))
  
  # Create list containing all output and return
  complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = method, "alpha" = alpha,
                                   "cluster.sizes" = cluster.sizes, "n.clusters" = n.clusters, "variance.parms" = var.parms, 
                                   "inputs" = difference, "model.estimates" = cps.model.est, "sim.data" = simulated.datasets), 
                              class = 'crtpwr')
  return(complete.output)
  }
