#' Power simulations for cluster-randomized trials: Stepped Wedge Design, Continuous Outcome.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for stepped wedge cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs power simulations for stepped wedge cluster-randomized controlled trials with continuous outcome.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatment arms, within-cluster variance, between-cluster variance, 
#' significance level, analytic method, progress updates, and simulated data 
#' set output may also be specified.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (equal cluster sizes) 
#' or a vector of length \code{nclusters} (user-defined size for each cluster) (required).
#' @param nclusters Number of clusters; accepts non-negative integer scalar (required).
#' @param difference Expected absolute difference in treatment effect between treatment and non-treatment groups;
#'  accepts numeric (required).
#' @param steps Number of crossover steps; accepts positive scalar (indicating the total number of steps; 
#' clusters per step is obtained by \code{nclusters / steps}) or a vector of non-negative integers corresponding 
#' either to the number of clusters to be crossed over at each time point (e.g c(2,4,4,2); nclusters = 10) or the 
#' cumulative number of clusters crossed over by a given time point (e.g. c(2,4,8,10); nclusters = 10) (required).
#' @param sigma Within-cluster variance; accepts non-negative numeric scalar (indicating equal within-cluster variances for both 
#' treatment groups) or a vector of length 2 specifying within-cluster variances for the non-treatment and treatment groups, 
#' respectively (required).
#' @param sigma_b Between-cluster variance; accepts non-negative numeric scalar (indicating equal 
#' between-cluster variances for both treatment groups) or a vector of length 2 specifying treatment-specific 
#' between-cluster variances (required).
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
#'   \item{variance.parms}{Data frame reporting ICC, within & between cluster variances for Treatment/Non-Treatment groups at each time point}
#'   \item{inputs}{Vector containing expected difference between groups based on user inputs}
#'   \item{means}{Data frame containing mean response values for each treatment group at each time point}
#'   \item{model.estimates}{Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)}
#'   \item{sim.data}{List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group),
#'                   "time.point" (Indicator for step; "t1" = time point 0) 
#'                   "clust" (Indicator for cluster), 
#'                   "period" (Indicator for at which step a cluster crosses over)}
#' }
#' 
#' @author Alexander R. Bogdan
#' 
#' @examples 
#' \dontrun{
#' normal.sw.rct = cps.sw.normal(nsim = 100, nsubjects = 50, nclusters = 30, 
#'                               difference = 3, steps = 5, sigma = 100, sigma_b = 30, 
#'                               alpha = 0.05, method = 'glmm', quiet = FALSE, 
#'                               all.sim.data = FALSE)
#' }
#'
#' @export

cps.sw.normal = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, difference = NULL, 
                         steps = NULL, sigma = NULL, sigma_b = NULL, alpha = 0.05, 
                         method = 'glmm', quiet = FALSE, all.sim.data = FALSE){
  
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
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS, DIFFERENCE
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
  if(length(nclusters) != 1){
    stop("NCLUSTERS must be a scalar (Total number of clusters)")
  }
  # Set sample sizes for each cluster (if not already specified)
  if(length(nsubjects) == 1){
    nsubjects[1:nclusters] = nsubjects
  }
  if(!length(nsubjects) %in% c(1, nclusters)){
    stop("NSUBJECTS must either be a scalar (indicating equal cluster sizes) or a vector of length NCLUSTERS (specifying cluster sizes for each cluster).")
  }
  
  # Validate STEPS
  if(!is.wholenumber(steps) || steps < 0){
    stop("All values supplied to STEPS must be non-negative integers")
  }
  if(length(steps) > 1){
    if(sum(steps) != nclusters && max(steps) != nclusters){
      stop("Total number of clusters specified by STEPS must either sum to NCLUSTERS or increase monotonically such that max(STEPS) == NCLUSTERS")
    }
  }
  if(length(steps) == 1){
    crossover.ind = 1:steps
    step.increment = nclusters / steps
    step.index = seq(from = step.increment, to = nclusters, by = step.increment)
    ### Need to figure out how to allocate clusters when nclusters %% steps != 0
  }
  # Create indexing vector for when SUM(STEPS) == NCLUSTERS & MAX(STEPS) == NCLUSTERS
  if(sum(steps) == nclusters){
    step.index = sapply(1:length(steps), function(x) sum(steps[0:x]))
    crossover.ind = 1:length(steps)
  }
  if(max(steps) == nclusters){
    crossover.ind = 1:length(steps)
    step.index = steps
  }
  
  # Create vector to store group means at each time point
  values.vector = cbind(c(rep(0, length(step.index) * 2)))

  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(difference) || difference < 0){
    stop("DIFFERENCE", min0.warning)
  }
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate SIGMA, SIGMA_B0, SIGMA_B1
  sigma_b.warning = " must be a scalar (equal between-cluster variance for both treatment and non-treatment groups) 
  or a vector of length 2, specifying unique between-cluster variances for the treatment and non-treatment groups."
  if(!is.numeric(sigma) || any(sigma < 0)){
    stop("All values supplied to SIGMA must be numeric values > 0")
  }
  if(!length(sigma) %in% c(1, 2)){
    stop("SIGMA must be a scalar (equal within-cluster variance for both treatment and non-treatment groups) 
         or a vector of length 2, specifying unique within-cluster variances for the treatment and non-treament groups.")
  }
  if(!is.numeric(sigma_b) || any(sigma_b < 0)){
    stop("All values supplied to SIGMA_B must be numeric values >= 0")
  }
  if(!length(sigma_b) %in% c(1, 2)){
    stop("SIGMA_B", sigma_b.warning)
  }
  # Set SIGMA & SIGMA_B (if not already set)
  if(length(sigma) == 1){
    sigma[2] = sigma
  }
  ### This (additive SIGMA_B) doesn't feel right or seem intuitive ###
  if(length(sigma_b) == 2){
    sigma_b[2] = sigma_b[1] + sigma_b[2]
  }
  if(length(sigma_b) == 1){
    sigma_b[2] = sigma_b
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
  
  # Create indicators for CLUSTER, STEP (period) & CROSSOVER (trt)
  clust = unlist(lapply(1:nclusters, function(x) rep(x, length.out = nsubjects[x])))
  period = NULL
  k = 1 # iterator
  for(i in 1:nclusters){
    if((i - 1) %in% step.index){
      k = k + 1
    }
    period = append(period, rep(crossover.ind[k], length.out = nsubjects[i]))
  }
  d = data.frame(period = period, clust = clust)
  
  # Create crossover indicators and melt output columns into single column
  for(j in 1:(length(crossover.ind) + 1)){
    d[[paste0("t", j)]] = ifelse(j > d[, 'period'], 1, 0)
  }
  sim.dat = tidyr::gather(d, key = 'time.point', value = 'trt', c(colnames(d)[-c(1,2)]))
  #sim.dat['time.point'] = as.numeric(gsub("t", "", sim.dat[, 'time.point']))
  sim.dat['y'] = 0
  
  ## Craete simulation & analysis loop
  for (i in 1:nsim){
    # Create vectors of cluster effects
    ntrt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b[1]))
    trt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b[2]))
    
    # Add subject specific effects & cluster effects
    for(j in 1:nclusters){
      # Assign non-treatment subject & cluster effects 
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 0, 
                            rnorm(sum(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 0), 0, sigma[1]) + 
                              ntrt.cluster.effects[j], 
                            sim.dat[, 'y'])
      # Assign treatment subject & cluster effects
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 1, 
                            rnorm(sum(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 1), difference, sigma[2]) + 
                              trt.cluster.effects[j], 
                            sim.dat[, 'y'])
    }
    # Add subject-specific error terms
    sim.dat['y'] = sim.dat['y'] + rnorm(nrow(sim.dat), 0, 1)
    
    # Calculate mean values for each group at each time point, for a given simulation
    iter.values = cbind(stats::aggregate(y ~ trt + time.point, data = sim.dat, mean)[, 3])
    values.vector = values.vector + iter.values
    
    # Store simulated data sets if ALL.SIM.DATA == TRUE 
    if(all.sim.data == TRUE){
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    ###################################################################
    ### DEV NOTE: Hussey & Hughes (2007) does not specify degrees of freedom for significance testing
    ###           GLMM DF = NCLUSTERS[Total # of clusters] - length(STEP.INDEX)[# of steps] - 2[# of treatment groups]
    ###################################################################
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + time.point + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = nclusters - length(step.index) - 2)
      est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt', 't value'])
      pval.vector = append(pval.vector, p.val)
    }
    
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y ~ trt + time.point, data = sim.dat,
                               id = clust, corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
    }
    
    # Update progress information
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
  
  ## Output objects
  # Create object containing summary statement
  summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, " Simulations: Stepped Wedge Design, Continuous Outcome")
  
  # Store simulation output in data frame
  cps.model.est = data.frame(Estimate = as.vector(unlist(est.vector)),
                             Std.err = as.vector(unlist(se.vector)),
                             Test.statistic = as.vector(unlist(stat.vector)),
                             p.value = as.vector(unlist(pval.vector)))
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  pval.power = sum(cps.model.est[, 'sig.val']) / nrow(cps.model.est)
  power.parms = data.frame(Power = round(pval.power, 3),
                           Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))

  # Create object containing treatment & time-specific differences
  values.vector = values.vector / nsim
  group.means = data.frame(Time.point = c(0, rep(1:(length(step.index) - 1), each = 2), length(step.index)), 
                           Treatment = c(0, rep(c(0, 1), length.out = (length(step.index) - 1) * 2), 1), 
                           Mean = round(values.vector, 3))
  
  # Create object containing cluster sizes
  cluster.sizes = nsubjects

  # Create object containing number of clusters
  n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), "Treatment" = c("n.clust" = nclusters[2])))

  # Create object containing variance parameters for each group at each time point
  var.parms = t(data.frame('Non.Treatment' = c('sigma' = sigma[1], 'sigma_b' = sigma_b[1]), 
                           'Treatment' = c('sigma' = sigma[2], 'sigma_b' = sigma_b[2])))

  # Create list containing all output (class 'crtpwr') and return
  complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = method, "alpha" = alpha,
                                   "cluster.sizes" = cluster.sizes, "n.clusters" = n.clusters, "variance.parms" = var.parms,
                                   "inputs" = difference, "means" = group.means, "model.estimates" = cps.model.est, "sim.data" = simulated.datasets),
                              class = 'crtpwr')

  return(complete.output)
}
