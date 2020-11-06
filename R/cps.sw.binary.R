#' Power simulations for cluster-randomized trials: Stepped Wedge Design, Binary Outcome
#'
#' This set of functions utilize iterative simulations to determine
#' approximate power for stepped wedge cluster-randomized controlled trials. Users
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#'
#' Runs power simulations for stepped wedge cluster-randomized controlled trials
#' with a binary outcome. The stepped wedge trial design is a type of cross-over
#' design in which clusters change treatments in waves. Initially all the
#' clusters recieve the same standard treatment, and at the end of the trial all
#' of the clusters will be recieving the treatment of interest. More than one
#' cluster can change treatments in a wave, but the order in which clusters
#' change treatments is randomly determined. The outcome of interest is assessed
#' in each cluster during each wave.
#'
#' Users must specify the desired number of simulations, number of subjects per
#' cluster, number of clusters per arm, expected absolute difference
#' between arms, within-cluster variance, between-cluster variance,
#' significance level, analytic method, progress updates, and simulated data
#' set output may also be specified.
#'
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (equal cluster sizes)
#' or a vector of length \code{nclusters} (user-defined size for each cluster) (required).
#' @param nclusters Number of clusters; accepts non-negative integer scalar (required).
#' @param p1 Expected probability of outcome in arm 1. Accepts scalar between 0 - 1 (required).
#' @param p2 Expected probability of outcome in arm 2. Accepts scalar between 0 - 1 (required).
#' @param steps Number of crossover steps; a baseline step (all clusters in arm 1) is assumed.
#' Accepts positive scalar (indicating the total number of steps; clusters per step is obtained by
#' \code{nclusters / steps}) or a vector of non-negative integers corresponding either to the number
#' of clusters to be crossed over at each time point (e.g c(2,4,4,2); nclusters = 10) or the cumulative
#' number of clusters crossed over by a given time point (e.g. c(2,4,8,10); nclusters = 10) (required).
#' @param sigma_b_sq Between-cluster variance; accepts non-negative numeric scalar (indicating equal
#' between-cluster variances for both arms) or a vector of length 2 specifying treatment-specific
#' between-cluster variances (required).
#' @param alpha Significance level. Default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; 
#' default is FALSE.
#' @param allSimData Option to output list of all simulated datasets; default = FALSE.
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge; default = FALSE.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' @param optmethod Option to fit with a different optimizer method (using the package \code{optimx}). 
#' Default is 'L-BFGS-B'.
#' @param seed Option to set.seed. Default is NULL.
#'
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations and simulation type
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power),
#'                "lower.95.ci" (Lower 95% confidence interval bound),
#'                "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting ICC, within & between cluster variances for both arms at each time point
#'   \item Vector containing expected difference between groups based on user inputs
#'   \item Data frame containing mean response values for each arm at each time point
#'   \item Matrix showing cluster crossover at each time point
#'   \item Data frame with columns:
#'                   "Estimate" (Estimate of treatment effect for a given simulation),
#'                   "Std.err" (Standard error for treatment effect estimate),
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)),
#'                   "p.value",
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing:
#'                   "y" (Simulated response value),
#'                   "trt" (Indicator for arm),
#'                   "time.point" (Indicator for step; "t1" = time point 0)
#'                   "clust" (Indicator for cluster),
#'                   "period" (Indicator for at which step a cluster crosses over)
#' }
#'  If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "clust" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#'
#' @examples
#'
#' # Estimate power for a trial with 3 steps and 12 clusters in arm 1
#' # (often the standard-of-care or 'control' arm) at the initiation of the study.
#' # Those clusters have 50 subjects each, with sigma_b_sq = 1.
#' # We have estimated arm outcome proportions of 0.1 and 0.2 in the first and second arms,
#' # respectively, and 100 simulated data sets analyzed by the GLMM method. Using seed = 123,
#' # the resulting power should be 0.8.
#'
#' \dontrun{
#' binary.sw.rct = cps.sw.binary(nsim = 100, nsubjects = 50, nclusters = 12,
#'                               p1 = 0.1, p2 = 0.2, steps = 3,
#'                               sigma_b_sq = 1, alpha = 0.05, method = 'glmm',
#'                               quiet = FALSE, allSimData = FALSE, seed = 123)
#' }
#'
#' @author Alexander R. Bogdan
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @export


cps.sw.binary = function(nsim = NULL,
                         nsubjects = NULL,
                         nclusters = NULL,
                         p1 = NULL,
                         p2 = NULL,
                         steps = NULL,
                         sigma_b_sq = NULL,
                         alpha = 0.05,
                         method = 'glmm',
                         quiet = FALSE,
                         allSimData = FALSE,
                         poorFitOverride = FALSE,
                         lowPowerOverride = FALSE,
                         timelimitOverride = TRUE,
                        optmethod= 'L-BFGS-B',
                         seed = NULL) {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }

  # Create vectors to collect iteration-specific values
  est.vector = vector("numeric", length = nsim)
  se.vector = vector("numeric", length = nsim)
  stat.vector = vector("numeric", length = nsim)
  pval.vector = vector("numeric", length = nsim)
  lmer.icc.vector = vector("numeric", length = nsim)
  converge = rep(NA, times=nsim)
  simulated.datasets = list()
  
  # Set start.time for progress iterator & initialize progress bar
  start.time = Sys.time()
  prog.bar =  progress::progress_bar$new(
    format = "(:spin) [:bar] :percent eta :eta",
    total = nsim,
    clear = FALSE,
    width = 100
  )
  prog.bar$tick(0)
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5)
    abs(x - round(x)) < tol
  
  # Define expit function
  expit = function(x)
    1 / (1 + exp(-x))
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS, DIFFERENCE
  if (length(steps) == 1) {
    if (!is.wholenumber(nclusters / steps)) {
      stop(
        "nclusters/steps must be a whole number. See documentation for steps parameter in '?clusterPower::cps.sw.binary'"
      )
    }
  }
  sim.data.arg.list = list(nsim, nclusters, nsubjects)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if (sum(sim.data.args) > 0) {
    stop(
      "NSIM, NCLUSTERS, NSUBJECTS & DIFFERENCE must all be specified. Please review your input values."
    )
  }
  min1.warning = " must be an integer greater than or equal to 1"
  if (!is.wholenumber(nsim) || nsim < 1) {
    stop(paste0("NSIM", min1.warning))
  }
  if (!is.wholenumber(nclusters) || nclusters < 1) {
    stop(paste0("NCLUSTERS", min1.warning))
  }
  if (!is.wholenumber(nsubjects) || nsubjects < 1) {
    stop(paste0("NSUBJECTS", min1.warning))
  }
  if (length(nclusters) != 1) {
    stop("NCLUSTERS must be a scalar (Total number of clusters)")
  }
  # Set sample sizes for each cluster (if not already specified)
  if (length(nsubjects) == 1) {
    nsubjects[1:nclusters] = nsubjects
  }
  if (!length(nsubjects) %in% c(1, nclusters)) {
    stop(
      "NSUBJECTS must either be a scalar (indicating equal cluster sizes) or a vector of length NCLUSTERS (specifying cluster sizes for each cluster)."
    )
  }
  
  # Validate P.NTRT & P.TRT
  min0.warning = " must be a numeric value between 0 - 1"
  if (p1 < 0 || p1 > 1) {
    stop("P.NTRT", min0.warning)
  }
  if (p2 < 0 || p2 > 1) {
    stop("P.TRT", min0.warning)
  }
  
  # Validate STEPS
  if (!is.wholenumber(steps) || steps < 0) {
    stop("All values supplied to STEPS must be non-negative integers")
  }
  if (length(steps) > 1) {
    if (sum(steps) != nclusters && max(steps) != nclusters) {
      stop(
        "Total number of clusters specified by STEPS must either sum to NCLUSTERS or increase monotonically such that max(STEPS) == NCLUSTERS"
      )
    }
  }
  if (length(steps) == 1) {
    crossover.ind = 1:steps
    step.increment = nclusters / steps
    step.index = seq(from = step.increment, to = nclusters, by = step.increment)
    ### Need to figure out how to allocate clusters when nclusters %% steps != 0
  }
  # Create indexing vector for when SUM(STEPS) == NCLUSTERS
  if (sum(steps) == nclusters) {
    step.index = sapply(1:length(steps), function(x)
      sum(steps[0:x]))
    crossover.ind = 1:length(steps)
  }
  if (max(steps) == nclusters) {
    crossover.ind = 1:length(steps)
    step.index = steps
  }
  
  # Create vector to store group means at each time point
  values.vector = cbind(c(rep(0, length(step.index) * 2)))
  
  # Validate sigma_b_sq
  sigma_b_sq.warning = " must be a scalar (equal between-cluster variance for both arms)
  or a vector of length 2, specifying unique between-cluster variances for both arms."
  if (!is.numeric(sigma_b_sq) || any(sigma_b_sq < 0)) {
    stop("All values supplied to sigma_b_sq must be numeric values >= 0")
  }
  if (!length(sigma_b_sq) %in% c(1, 2)) {
    stop("sigma_b_sq", sigma_b_sq.warning)
  }
  # Set sigma_b_sq (if not already set)
  # Note: If user-defined, sigma_b_sq[2] is additive
  if (length(sigma_b_sq) == 2) {
    sigma_b_sq[2] = sigma_b_sq[1] + sigma_b_sq[2]
  }
  if (length(sigma_b_sq) == 1) {
    sigma_b_sq[2] = sigma_b_sq
  }
  
  # Validate ALPHA, METHOD, QUIET, allSimData
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA", min0.warning)
  }
  if (!is.element(method, c('glmm', 'gee'))) {
    stop(
      "METHOD must be either 'glmm' (Generalized Linear Mixed Model)
         or 'gee'(Generalized Estimating Equation)"
    )
  }
  if (!is.logical(quiet)) {
    stop(
      "QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)"
    )
  }
  if (!is.logical(allSimData)) {
    stop(
      "allSimData must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output"
    )
  }
  browser()
  # Calculate ICC1 (sigma_b_sq / (sigma_b_sq + pi^2/3))
  icc1 = mean(sapply(1:2, function(x)
    sigma_b_sq[x] / (sigma_b_sq[x] + pi ^ 2 / 3)))
  
  # Create indicators for CLUSTER, STEP (period) & CROSSOVER (trt)
  clust = unlist(lapply(1:nclusters, function(x)
    rep(x, length.out = nsubjects[x])))
  period = NULL
  k = 1 # iterator
  for (i in 1:nclusters) {
    if ((i - 1) %in% step.index) {
      k = k + 1
    }
    period = append(period, rep(crossover.ind[k], length.out = nsubjects[i]))
  }
  d = data.frame(period = period, clust = clust)
  
  # Create crossover indicators and melt output columns into single column
  for (j in 1:(length(crossover.ind) + 1)) {
    d[[paste0("t", j)]] = ifelse(j > d[, 'period'], 1, 0)
  }
  sim.dat = tidyr::gather(d, key = 'time.point', value = 'trt', c(colnames(d)[-c(1, 2)]))
  sim.dat['y'] = 0
  
  # Calculate log odds for each group
  logit.p1 = log(p1 / (1 - p1))
  logit.p2 = log(p2 / (1 - p2))
  
  ## Create simulation & analysis loop
  for (i in 1:nsim) {
    # Create vectors of cluster effects
    ntrt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b_sq[1]))
    trt.cluster.effects = stats::rnorm(nclusters, mean = 0, sd = sqrt(sigma_b_sq[2]))
    
    # Add subject specific effects & cluster effects
    for (j in 1:nclusters) {
      # Assign arm 1 subject & cluster effects
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == j &
                              sim.dat[, 'trt'] == 0,
                            stats::rbinom(
                              sum(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 0),
                              1,
                              expit(logit.p1 +
                                      ntrt.cluster.effects[j] +
                                      stats::rnorm(sum(
                                        sim.dat[, 'clust'] == j &
                                          sim.dat[, 'trt'] == 0
                                      )))
                            ),
                            sim.dat[, 'y'])
      # Assign arm 2 subject & cluster effects
      sim.dat['y'] = ifelse(sim.dat[, 'clust'] == j &
                              sim.dat[, 'trt'] == 1,
                            stats::rbinom(
                              sum(sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 1),
                              1,
                              expit(logit.p2 +
                                      trt.cluster.effects[j] +
                                      stats::rnorm(sum(
                                        sim.dat[, 'clust'] == j & sim.dat[, 'trt'] == 1
                                      )))
                            ),
                            sim.dat[, 'y'])
    }
    
    # Calculate LMER.ICC (lmer: sigma_b_sq / (sigma_b_sq + sigma))
    lmer.mod = lme4::lmer(y ~ trt + time.point + (1 |
                                                    clust), data = sim.dat)
    lmer.vcov = as.data.frame(lme4::VarCorr(lmer.mod))[, 4]
    lmer.icc.vector =  append(lmer.icc.vector, lmer.vcov[1] / (lmer.vcov[1] + lmer.vcov[2]))
    
    # Calculate mean values for each group at each time point, for a given simulation
    iter.values = cbind(stats::aggregate(y ~ trt + time.point, data = sim.dat, mean)[, 3])
    values.vector = values.vector + iter.values
    
    # Store simulated data sets if allSimData == TRUE
    if (allSimData == TRUE) {
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    sim.dat$period <- as.factor(sim.dat$period)
    sim.dat$clust <- as.factor(sim.dat$clust)
    sim.dat$time.point <- as.factor(sim.dat$time.point)
    sim.dat$trt <- as.factor(sim.dat$trt)
    sim.dat$y <- as.factor(sim.dat$y)
    
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      # Option to use optimizerSearch
      if (i == 1 & optmethod == "auto") {
        my.mod = lme4::glmer(
          y ~ trt + time.point + (1 | clust),
          data = sim.dat,
          family = stats::binomial(link = 'logit')
          )
        optmethod <- optimizerSearch(my.mod)
      } 
      my.mod = lme4::glmer(
        y ~ trt + time.point + (1 | clust),
        data = sim.dat,
        family = stats::binomial(link = 'logit'),
        control = lme4::glmerControl(
          optimizer = "optimx",
          calc.derivs = TRUE,
          optCtrl = list(
            method = optmethod,
            starttests = FALSE,
            kkt = FALSE
          )
        )
      )
      glmm.values = summary(my.mod)$coefficient
      est.vector[i] = glmm.values['trt', 'Estimate']
      se.vector[i] = glmm.values['trt', 'Std. Error']
      stat.vector[i] = glmm.values['trt', 'z value']
      pval.vector[i] = glmm.values['trt', 'Pr(>|z|)']
      converge[i] = is.null(my.mod@optinfo$conv$lme4$messages)
    }
    
    # Fit GEE (geeglm)
    if (method == 'gee') {
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(
        y ~ trt + time.point,
        data = sim.dat,
        family = stats::binomial(link = 'logit'),
        id = clust,
        corstr = "exchangeable"
      )
      gee.values = summary(my.mod)$coefficients
      est.vector[i] = gee.values['trt', 'Estimate']
      se.vector[i] = gee.values['trt', 'Std.err']
      stat.vector[i] = gee.values['trt', 'Wald']
      pval.vector[i] = gee.values['trt', 'Pr(>|W|)']
      converge[i] = ifelse(summary(my.mod)$error == 0, TRUE, FALSE)
    }
    
    # option to stop the function early if fits are singular
    if (poorFitOverride == FALSE &&
        converge[i] == FALSE && i > 50) {
      if (sum(converge == FALSE, na.rm = TRUE) > (nsim * .25)) {
        stop("more than 25% of simulations are singular fit: check model specifications")
      }
    }
    
    # stop the loop if power is <0.5
    if (lowPowerOverride == FALSE && i > 50 && (i %% 10 == 0)) {
      sig.val.temp <-
        ifelse(pval.vector < alpha, 1, 0)
      pval.power.temp <- sum(sig.val.temp, na.rm = TRUE) / i
      if (pval.power.temp < 0.5) {
        stop(
          paste(
            "Calculated power is < ",
            pval.power.temp,
            ". Set lowPowerOverride == TRUE to run the simulations anyway.",
            sep = ""
          )
        )
      }
    }
    
    # Update progress information
    if (i == 1) {
      avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      time.est = avg.iter.time * (nsim - 1) / 60
      hr.est = time.est %/% 60
      min.est = round(time.est %% 60, 0)
      if (min.est > 2 && timelimitOverride == FALSE) {
        stop(paste0(
          "Estimated completion time: ",
          hr.est,
          'Hr:',
          min.est,
          'Min'
        ))
      }
      if (quiet == FALSE) {
        message(
          paste0(
            'Begin simulations :: Start Time: ',
            Sys.time(),
            ' :: Estimated completion time: ',
            hr.est,
            'Hr:',
            min.est,
            'Min'
          )
        )
      }
    }
    # Iterate progress bar
    prog.bar$update(i / nsim)
    Sys.sleep(1 / 100)
    
    if (i == nsim) {
      total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      hr.est = total.est %/% 3600
      min.est = total.est %/% 60
      sec.est = round(total.est %% 60, 0)
      message(
        paste0(
          "Simulations Complete! Time Completed: ",
          Sys.time(),
          "\nTotal Runtime: ",
          hr.est,
          'Hr:',
          min.est,
          'Min:',
          sec.est,
          'Sec'
        )
      )
    }
  }
  
  ## Output objects
  # Create object containing summary statement
  summary.message = paste0(
    "Monte Carlo Power Estimation based on ",
    nsim,
    " Simulations: Stepped Wedge Design, Binary Outcome"
  )
  
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model',
                       gee = 'Generalized Estimating Equation')
  
  # Store simulation output in data frame
  cps.model.est = data.frame(
    Estimate = as.vector(unlist(est.vector)),
    Std.err = as.vector(unlist(se.vector)),
    Test.statistic = as.vector(unlist(stat.vector)),
    p.value = as.vector(unlist(pval.vector)),
    converge = as.vector(converge)
  )
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  cps.model.temp <- dplyr::filter(cps.model.est, converge == TRUE)
  power.parms <- confintCalc(alpha = alpha,
                             p.val = cps.model.temp[, 'p.value'])
  
  # Create object containing treatment & time-specific differences
  values.vector = values.vector / nsim
  group.means = data.frame(
    Time.point = c(0, rep(1:(
      length(step.index) - 1
    ), each = 2), length(step.index)),
    Arm = c(0, rep(c(0, 1), length.out = (
      length(step.index) - 1
    ) * 2), 1),
    Values = round(values.vector, 3)
  )
  
  # Create object containing expected arm 1 and arm 2 probabilities
  group.probs = data.frame("Outcome.Probabilities" = c("Arm.1" = p1, "Arm.2" = p2))
  
  # Create object containing cluster sizes
  cluster.sizes = nsubjects
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm.1" = c("n.clust" = nclusters[1]),
    "Arm.2" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing estimated ICC values
  ICC = round(t(data.frame(
    'P_h' = c('ICC' = icc1),
    'lmer' = c('ICC' = mean(lmer.icc.vector))
  )), 3)
  
  # Create object containing variance parameters for each group at each time point
  var.parms = t(data.frame(
    'Arm.1' = c('sigma_b_sq' = sigma_b_sq[1]),
    'Arm.2' = c('sigma_b_sq' = sigma_b_sq[2])
  ))
  
  # Create crossover matrix output object
  crossover.mat = apply(as.matrix(c(0, step.index)), 1,
                        function(x)
                          c(rep(1, length.out = x), rep(0, length.out = nclusters - x)))
  
  # Create list containing all output (class 'crtpwr') and return
  complete.output = structure(
    list(
      "overview" = summary.message,
      "nsim" = nsim,
      "power" = power.parms,
      "method" = long.method,
      "alpha" = alpha,
      "cluster.sizes" = cluster.sizes,
      "n.clusters" = n.clusters,
      "variance.parms" = var.parms,
      "inputs" = group.probs,
      "means" = group.means,
      "icc" = ICC,
      "model.estimates" = cps.model.est,
      "sim.data" = simulated.datasets,
      "crossover.matrix" = crossover.mat
    ),
    class = 'crtpwr'
  )
  
  return(complete.output)
}
