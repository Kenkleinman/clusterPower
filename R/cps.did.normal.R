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
#' @param nclusters Number of clusters per group; accepts integer scalar or vector of length 2 for unequal number 
#' of clusters per treatment group (required)
#' @param difference Expected absolute difference in treatment effect between time points; accepts numeric (required).
#' @param sigma Within-cluster variance; accepts numeric scalar (indicating equal within-cluster variances for both 
#' treatment groups at both time points) or vector of length 4 specifying within-cluster variance for each treatment 
#' group at each time point.
#' @param sigma_b0 Pre-treatment (time == 0) between-cluster variance; accepts numeric scalar (indicating equal 
#' between-cluster variances for both treatment groups) or a vector of length 2 specifying treatment-specific 
#' between-cluster variances
#' @param sigma_b1 Post-treatment (time == 1) between-cluster variance; accepts numeric scalar (indicating equal 
#' between-cluster variances for both treatment groups) or a vector of length 2 specifying treatment-specific 
#' between-cluster variances. For data simulation, SIGMA_B1 is added to SIGMA_B0, such that if SIGMA_B0 = 5 
#' and SIGMA_B1 = 2, the between-cluster variance at time == 1 equals 7. Default = 0.
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
#'   \item{differences}{Data frame with columns: 
#'                   "Period" (Pre/Post-treatment indicator), 
#'                   "Treatment" (Treatment group indicator), 
#'                   "Value" (Mean response value)}
#'   \item{model.estimates}{Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)}
#'   \item{sim.data}{List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster), 
#'                   "period" (Indicator for time point)}
#' }
#' 
#' @author Alexander R. Bogdan
#' 
#' @examples 
#' \dontrun{
#' normal.did.rct = cps.did.normal(nsim = 100, nsubjects = 50, nclusters = 30, 
#'                                 difference = 3, ICC = 0.2, sigma = 100, alpha = 0.05, 
#'                                 method = 'glmm', quiet = FALSE, all.sim.data = FALSE)
#' }
#'
#' @export

cps.did.normal = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, difference = NULL, 
                          sigma = NULL, sigma_b0 = NULL, sigma_b1 = 0, alpha = 0.05, 
                          method = 'glmm', quiet = FALSE, all.sim.data = FALSE){
  
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  values.vector = cbind(c(0, 0, 0, 0))
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
  
  # Validate SIGMA, SIGMA_B0, SIGMA_B1
  sigma_b.warning = " must be a scalar (equal between-cluster variance for both treatment groups) or a vector of length 2, 
         specifying between-cluster variances for each treatment group"
  if(!is.numeric(sigma) || any(sigma < 0)){
    stop("All values supplied to SIGMA must be numeric values > 0")
  }
  if(!length(sigma) %in% c(1,4)){
    stop("SIGMA must be a scalar (equal within-cluster variance for both treatment groups at both time points) 
         or a vector of length 4, specifying within-cluster variances for each treatment group at each time point")
  }
  if(!is.numeric(sigma_b0) || any(sigma_b0 < 0)){
    stop("All values supplied to SIGMA_B0 must be numeric values > 0")
  }
  if(!length(sigma_b0) %in% c(1,2)){
    stop("SIGMA_B0", sigma_b.warning)
  }
  if(!length(sigma_b1) %in% c(1,2)){
    stop("SIGMA_B1", sigma_b.warning)
  }
  if(!is.numeric(sigma_b1) || any(sigma_b1 < 0)){
    stop("All values supplied to SIGMA_B1 must be numeric values >= 0")
  }
  # Set SIGMA, SIGMA_B0 & SIGMA_B1 (if not already set)
  if(length(sigma) == 1){
    sigma = rep(sigma, 4)
  }
  if(length(sigma_b0) == 1){
    sigma_b0[2] = sigma_b0
  }
  if(length(sigma_b1) == 1){
    sigma_b1[2] = sigma_b1
  }
  sigma_b1 = sigma_b1 + sigma_b0
  
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
  
  # Create indicators for PERIOD, TRT & CLUST
  period = rep(0:1, each = sum(nsubjects))
  trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])), 
          rep(1, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))

  # Create simulation loop
  for(i in 1:nsim){
  ### Generate simulated data
    ## TIME == 0
    # Generate between-cluster effects for non-treatment and treatment
    randint.ntrt.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b0[1]))
    randint.trt.0 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b0[2]))

    # Create non-treatment y-value
    y0.ntrt.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.ntrt.0[x], length.out = nsubjects[x])))
    y0.ntrt.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[1]))))
    y0.ntrt.pre = y0.ntrt.bclust + y0.ntrt.wclust + stats::rnorm(nsubjects[1:nclusters[1]])

    # Create treatment y-value
    y0.trt.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.trt.0[x], length.out = nsubjects[nclusters[1] + x])))
    y0.trt.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                              function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[2]))))
    y0.trt.pre = y0.trt.bclust + y0.trt.wclust + stats::rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])

    ## TIME == 1
    # Generate between-cluster effects for non-treatment and treatment
    randint.ntrt.1 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b1[1]))
    randint.trt.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b1[2]))

    # Create non-treatment y-value
    y1.ntrt.bclust = unlist(lapply(1:nclusters[1], function(x) rep(randint.ntrt.1[x], length.out = nsubjects[x])))
    y1.ntrt.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x) stats::rnorm(x, mean = 0, sd = sqrt(sigma[3]))))
    y1.ntrt.post = y1.ntrt.bclust + y1.ntrt.wclust + stats::rnorm(nsubjects[1:nclusters[1]])

    # Create treatment y-value
    y1.trt.bclust = unlist(lapply(1:nclusters[2], function(x) rep(randint.trt.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.trt.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                              function(x) stats::rnorm(x, mean = difference, sd = sqrt(sigma[4]))))
    y1.trt.post = y1.trt.bclust + y1.trt.wclust + stats::rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])

    # Create single response vector
    y = c(y0.ntrt.pre, y0.trt.pre, y1.ntrt.post, y1.trt.post)

    # Create data frame for simulated dataset
    sim.dat = data.frame(y = y, trt = trt, clust = clust, period = period)
    if(all.sim.data == TRUE){
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # Calculate mean values for given simulation
    iter.values = cbind(stats::aggregate(y ~ trt + period, data = sim.dat, mean)[, 3])
    values.vector = values.vector + iter.values
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + period + trt:period + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt:period', 't value']), df = sum(nclusters) - 2)
      est.vector = append(est.vector, glmm.values['trt:period', 'Estimate'])
      se.vector = append(se.vector, glmm.values['trt:period', 'Std. Error'])
      stat.vector = append(stat.vector, glmm.values['trt:period', 't value'])
      pval.vector = append(pval.vector, p.val)
    }
    
    # Fit GEE (geeglm)
    if(method == 'gee'){
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y ~ trt + period + trt:period, data = sim.dat,
                               id = clust, corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt:period', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt:period', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt:period', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt:period', 'Pr(>|W|)'])
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
                           Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  
  # Create object containing treatment & time-specific differences
  values.vector = values.vector / nsim
  differences = data.frame(Period = c(0,0,1,1), Treatment = c(0,1,0,1), Values = round(values.vector, 3))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Non.Treatment' = nsubjects[1:nclusters[1]], 
                       'Treatment' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), "Treatment" = c("n.clust" = nclusters[2])))
  
  # Create object containing variance parameters for each group at each time point
  var.parms = list("Time.Point.0" = data.frame('Non.Treatment' = c("sigma" = sigma[1], "sigma_b" = sigma_b0[1]), 
                                                'Treatment' = c("sigma" = sigma[2], "sigma_b" = sigma_b0[2])), 
                   "Time.Point.1" = data.frame('Non.Treatment' = c("sigma" = sigma[3], "sigma_b" = sigma_b1[1]), 
                                            'Treatment' = c("sigma" = sigma[4], "sigma_b" = sigma_b1[2])))
  
  # Create list containing all output (class 'crtpwr') and return
  complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = method, "alpha" = alpha,
                                   "cluster.sizes" = cluster.sizes, "n.clusters" = n.clusters, "variance.parms" = var.parms, 
                                   "inputs" = difference, "model.estimates" = cps.model.est, "sim.data" = simulated.datasets, 
                                   "differences" = differences),
                              class = 'crtpwr')
  return(complete.output)
  }
