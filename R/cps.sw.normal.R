###############################################################################
### PROJECT: cps.sw.normal
### AUTHOR: Alexander R. Bogdan
### DATE CREATED: 04/30/18
### DATE LAST MODIFIED: 05/03/18
### DEV NOTES: 
## 1) The below code is in ALPHA. The simulation engine and analytic engines are FUNCTIONAL.
## 2) Revisit validation for STEPS, SIGMA_B0, SIGMA_B1
## 3) Figure out how to incorporate user-specified cluster corssover at a given time point
## 4) Figure out how to distribute clusters when NCLUSTER %% STEPS != 0
## 5) Decide how to change output and reconcile this with SUMMARY & PRINT methods
###############################################################################


cps.sw.normal = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, difference = NULL, 
                         steps = NULL, sigma = NULL, sigma_b0 = NULL, sigma_b1 = 0, alpha = 0.05, 
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
  if(!is.wholenumber(nsubjects) || steps < 1){
    stop(paste0("STEPS", min1.warning))
  }
  if(length(steps) == 1){
    steps = 1:steps
  }
  step.increment = nclusters / steps
  ################################################
  ### NEED TO ADD MORE VALIDATION HERE
  ################################################
  
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
  if(!is.numeric(sigma_b0) || any(sigma_b0 < 0)){
    stop("All values supplied to SIGMA_B0 must be numeric values > 0")
  }
  if(!length(sigma_b0) %in% c(1, 2)){
    stop("SIGMA_B0", sigma_b.warning)
  }
  if(!length(sigma_b1) %in% c(1, 2)){
    stop("SIGMA_B1", sigma_b.warning)
  }
  if(!is.numeric(sigma_b1) || any(sigma_b1 < 0)){
    stop("All values supplied to SIGMA_B1 must be numeric values >= 0")
  }
  # Set SIGMA, SIGMA_B0 & SIGMA_B1 (if not already set)
  if(length(sigma) == 1){
    sigma[2] = sigma
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
  
  # Create indicators for CLUSTER, STEP (period) & CROSSOVER (trt)
  clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))
  period = NULL
  k = 0
  for(i in 1:nclusters){
    if(i %% step.inc != 0 & (i - 1) %% step.inc == 0){
      k = k + 1
    }
    period = append(period, rep(steps[k], length.out = nsubjects[i]))
  }
  d = data.frame(period = period, clust = clust)
  
  # Create crossover indicators and melt output columns into single column
  for(j in 1:(length(steps) + 1)){
    d[[paste0("t", j)]] = ifelse(j > d[, 'period'], 1, 0)
  }
  sim.dat = tidyr::gather(d, key = 'time.point', value = 'trt', c(colnames(d)[-c(1,2)]))
  sim.dat['y'] = 0
  
  ### Loop
  for (i in 1:nsim){
    # Create vectors of cluster effects
    ntrt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b[1]))
    trt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b[2]))
    
    # Add subject specific effects & cluster effects
    for(i in 1:nclusters){
      # Assign non-treatment effects
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 0, 
                            rnorm(sum(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 0), 0, sigma[1]) + 
                              ntrt.clust.effect[i] + 
                              rnorm(sum(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 0), 0, 1), 
                            sim.dat[, 'y'])
      # Assign treatment effects
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 1, 
                            rnorm(sum(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 1), trt.effect, sigma[2]) + 
                              trt.clust.effect[i] + 
                              rnorm(sum(sim.dat[, 'clust'] == i & sim.dat[, 'trt'] == 1), 0, 1), 
                            sim.dat[, 'y'])
    }
    
    # Store simulated data sets if ALL.SIM.DATA == TRUE 
    if(all.sim.data == TRUE){
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # Calculate mean values for given simulation
    #iter.values = cbind(stats::aggregate(y ~ trt + period, data = sim.dat, mean)[, 3])
    #values.vector = values.vector + iter.values
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + time.point + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = sum(nclusters) - 2) # Degrees of freedom for p-value?
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
  
  ##################################
  ### BEGIN WORK FROM HERE DOWN
  ### Need to reconcile output with SUMMARY method
  ##################################
  
  
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