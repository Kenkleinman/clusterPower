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
#' cluster, number of clusters per arm, expected absolute difference
#' between arms, two of the following: ICC, within-cluster variance, or
#' between-cluster variance; significance level, analytic method, progress updates,
#' and simulated data set output may also be specified.
#'
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per arm; accepts either a scalar (equal cluster sizes, both groups),
#' a vector of length two (equal cluster sizes within groups), or a vector of length \code{sum(nclusters)}
#' (unequal cluster sizes within groups) (required).
#' @param nclusters Number of clusters per group; accepts integer scalar or vector of length 2 for unequal number
#' of clusters per arm (required)
#' @param mu Expected mean of arm 1; accepts numeric (required).
#' @param mu2 Expected mean of arm 2; accepts numeric (required).
#' @param sigma_sq Within-cluster variance; accepts numeric scalar (indicating equal within-cluster variances for both
#' arms at both time points) or vector of length 4 specifying within-cluster variance for each arm at each time point.
#' @param sigma_b_sq0 Pre-treatment (time == 0) between-cluster variance; accepts numeric scalar (indicating equal
#' between-cluster variances for both arms) or a vector of length 2 specifying treatment-specific
#' between-cluster variances
#' @param sigma_b_sq1 Post-treatment (time == 1) between-cluster variance; accepts numeric scalar (indicating equal
#' between-cluster variances for both arm) or a vector of length 2 specifying treatment-specific
#' between-cluster variances. For data simulation, sigma_b_sq1 is added to sigma_b_sq0, such that if sigma_b_sq0 = 5
#' and sigma_b_sq1 = 2, the between-cluster variance at time == 1 equals 7. Default = 0.
#' @param alpha Significance level. Default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge; default = FALSE.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param allSimData Option to output list of all simulated datasets; default = FALSE.
#' @param nofit Option to skip model fitting and analysis and only return the simulated data.
#' Default = \code{FALSE}.
#' @param seed Option to set the seed. Default is NA.
#'
#' @return A list with the following components:
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
#'   \item Data frame reporting ICC, within & between cluster variances
#'   for both arms at each time point
#'   \item Vector containing expected difference between groups based on user inputs
#'   \item Data frame with columns:
#'                   "Period" (Pre/Post-treatment indicator),
#'                   "Arm.2" (arm indicator),
#'                   "Value" (Mean response value)
#'   \item Data frame with columns:
#'                   "Estimate" (Estimate of treatment effect for a given simulation),
#'                   "Std.err" (Standard error for treatment effect estimate),
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)),
#'                   "p.value",
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing:
#'                   "y" (Simulated response value),
#'                   "trt" (Indicator for arm),
#'                   "clust" (Indicator for cluster),
#'                   "period" (Indicator for time point)
#' }
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#'
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#'
#' @examples
#'
#' # Estimate power for a trial with 6 clusters in arm 1 and 6 clusters in arm 2,
#' # those clusters having 120 subjects each, with sigma_sq = 1. Estimated
#' # arm means are 1 and 0.48 in the first and second arms, respectively, and we use
#' # 100 simulated data sets analyzed by the GLMM method. The resulting estimated
#' # power (for seed = 123) should be 0.81.
#'
#' \dontrun{
#' normal.did.rct = cps.did.normal(nsim = 100, nsubjects = 120, nclusters = 6,
#'                                 mu = 1, mu2 = 0.48, sigma_sq = 1, alpha = 0.05,
#'                                 sigma_b_sq0 = 0.1, method = 'glmm', seed = 123)
#' }
#'
#' @author Alexander R. Bogdan
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @export



cps.did.normal = function(nsim = NULL,
                          nsubjects = NULL,
                          nclusters = NULL,
                          mu = 0,
                          mu2 = NULL,
                          sigma_sq = NULL,
                          sigma_b_sq0 = NULL,
                          sigma_b_sq1 = 0,
                          alpha = 0.05,
                          method = 'glmm',
                          poorFitOverride = FALSE,
                          lowPowerOverride = FALSE,
                          timelimitOverride = TRUE,
                          quiet = FALSE,
                          allSimData = FALSE,
                          seed = NA,
                          nofit = FALSE) {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create vectors to collect iteration-specific values
  est.vector = vector("numeric", length = nsim)
  se.vector = vector("numeric", length = nsim)
  stat.vector = vector("numeric", length = nsim)
  pval.vector = vector("numeric", length = nsim)
  converge = vector(mode = "logical", length = nsim)
  values.vector = cbind(c(0, 0, 0, 0))
  simulated.datasets = list()
  
  # initialize progress bar
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
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS, mu and mu2
  sim.data.arg.list = list(nsim, nclusters, nsubjects, mu, mu2)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if (sum(sim.data.args) > 0) {
    stop(
      "nsim, nclusters, nsubjects, mu, and mu2 must all be specified. Please review your input values."
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
  if (length(nclusters) > 2) {
    stop(
      "NCLUSTERS can only be a vector of length 1 (equal # of clusters per group) or 2 (unequal # of clusters per group)"
    )
  }
  # Set cluster sizes for arm 2 arm (if not already specified)
  if (length(nclusters) == 1) {
    nclusters[2] = nclusters[1]
  }
  # Set sample sizes for each cluster (if not already specified)
  if (length(nsubjects) == 1) {
    nsubjects[1:sum(nclusters)] = nsubjects
  }
  if (length(nsubjects) == 2) {
    nsubjects = c(rep(nsubjects[1], nclusters[1]), rep(nsubjects[2], nclusters[2]))
  }
  if (nclusters[1] == nclusters[2] &&
      length(nsubjects) == nclusters[1]) {
    nsubjects = rep(nsubjects, 2)
  }
  if (length(nclusters) == 2 &&
      length(nsubjects) != 1 &&
      length(nsubjects) != sum(nclusters)) {
    stop(
      "A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS"
    )
  }
  
  # Validate mu and mu2, ALPHA
  min0.warning = "must be numeric values"
  if (!is.numeric(mu) || !is.numeric(mu2)) {
    stop("mu & mu2", min0.warning)
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate sigma_sq, sigma_b_sq0, sigma_b_sq1
  sigma_b_sq.warning = " must be a scalar (equal between-cluster variance for both arms) or a vector of length 2,
         specifying between-cluster variances for each arm"
  if (!is.numeric(sigma_sq) || any(sigma_sq < 0)) {
    stop("All values supplied to sigma_sq must be numeric values > 0")
  }
  if (!length(sigma_sq) %in% c(1, 4)) {
    stop(
      "sigma_sq must be a scalar (equal within-cluster variance for both arms at both time points)
         or a vector of length 4, specifying within-cluster variances for each arm at each time point"
    )
  }
  if (!is.numeric(sigma_b_sq0) || any(sigma_b_sq0 < 0)) {
    stop("All values supplied to sigma_b_sq0 must be numeric values > 0")
  }
  if (!length(sigma_b_sq0) %in% c(1, 2)) {
    stop("sigma_b_sq0", sigma_b_sq.warning)
  }
  if (!length(sigma_b_sq1) %in% c(1, 2)) {
    stop("sigma_b_sq1", sigma_b_sq.warning)
  }
  if (!is.numeric(sigma_b_sq1) || any(sigma_b_sq1 < 0)) {
    stop("All values supplied to sigma_b_sq1 must be numeric values >= 0")
  }
  # Set sigma_sq, sigma_b_sq0 & sigma_b_sq1 (if not already set)
  if (length(sigma_sq) == 1) {
    sigma_sq = rep(sigma_sq, 4)
  }
  if (length(sigma_b_sq0) == 1) {
    sigma_b_sq0[2] = sigma_b_sq0
  }
  if (length(sigma_b_sq1) == 1) {
    sigma_b_sq1[2] = sigma_b_sq1
  }
  sigma_b_sq1 = sigma_b_sq1 + sigma_b_sq0
  
  # Validate METHOD, QUIET, allSimData
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
  
  # Create indicators for PERIOD, TRT & CLUST
  period = rep(0:1, each = sum(nsubjects))
  trt = c(rep(1, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(2, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  
  start.time = Sys.time()
  
  # Create simulation loop
  for (i in 1:nsim) {
    ### Generate simulated data
    ## TIME == 0
    # Generate between-cluster effects for arm 1 and arm 2
    randint.ntrt.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq0[1]))
    randint.trt.0 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq0[2]))
    
    # Create arm 1 y-value
    y0.ntrt.bclust = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.ntrt.0[x], length.out = nsubjects[x])))
    y0.ntrt.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x)
      stats::rnorm(
        x, mean = 0, sd = sqrt(sigma_sq[1])
      )))
    y0.ntrt.pre = y0.ntrt.bclust + y0.ntrt.wclust + stats::rnorm(nsubjects[1:nclusters[1]])
    
    # Create arm 2 y-value
    y0.trt.bclust = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.trt.0[x], length.out = nsubjects[nclusters[1] + x])))
    y0.trt.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                                  function(x)
                                    stats::rnorm(
                                      x, mean = 0, sd = sqrt(sigma_sq[2])
                                    )))
    y0.trt.pre = y0.trt.bclust + y0.trt.wclust + stats::rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
    
    ## TIME == 1
    # Generate between-cluster effects for arm 1 and arm 2
    randint.ntrt.1 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq1[1]))
    randint.trt.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq1[2]))
    
    # Create arm 1 y-value
    y1.ntrt.bclust = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.ntrt.1[x], length.out = nsubjects[x])))
    y1.ntrt.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x)
      stats::rnorm(
        x, mean = mu, sd = sqrt(sigma_sq[3])
      )))
    y1.ntrt.post = y1.ntrt.bclust + y1.ntrt.wclust + stats::rnorm(nsubjects[1:nclusters[1]])
    
    # Create arm 2 y-value
    y1.trt.bclust = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.trt.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.trt.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                                  function(x)
                                    stats::rnorm(
                                      x, mean = mu2, sd = sqrt(sigma_sq[4])
                                    )))
    y1.trt.post = y1.trt.bclust + y1.trt.wclust + stats::rnorm(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
    
    # Create single response vector
    y = c(y0.ntrt.pre, y0.trt.pre, y1.ntrt.post, y1.trt.post)
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(
      y = y,
      trt = trt,
      clust = clust,
      period = period
    )
    if (allSimData == TRUE) {
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # option to return simulated data only
    if (nofit == TRUE) {
      if (!exists("nofitop")) {
        nofitop <- data.frame(
          period = sim.dat['period'],
          cluster = sim.dat['clust'],
          arm = sim.dat['trt'],
          y1 = sim.dat["y"]
        )
      } else {
        nofitop[, length(nofitop) + 1] <- sim.dat["y"]
      }
      if (length(nofitop) == (nsim + 3)) {
        temp1 <- seq(1:nsim)
        temp2 <- paste0("y", temp1)
        colnames(nofitop) <- c('period', 'cluster', 'arm', temp2)
      }
      if (length(nofitop) != (nsim + 3)) {
        next()
      }
      return(nofitop)
    }
    
    # Calculate mean values for given simulation
    iter.values = cbind(stats::aggregate(y ~ trt + period, data = sim.dat, mean)[, 3])
    values.vector = values.vector + iter.values
    
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      my.mod = lme4::lmer(y ~ trt + period + trt:period + (1 |
                                                             clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt:period', 't value']), df = sum(nclusters) - 2)
      est.vector[i] = glmm.values['trt:period', 'Estimate']
      se.vector[i] = glmm.values['trt:period', 'Std. Error']
      stat.vector[i] = glmm.values['trt:period', 't value']
      pval.vector[i] = p.val
      converge[i] = is.null(my.mod@optinfo$conv$lme4$messages)
    }
    
    # Fit GEE (geeglm)
    if (method == 'gee') {
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(
        y ~ trt + period + trt:period,
        data = sim.dat,
        id = clust,
        corstr = "exchangeable"
      )
      gee.values = summary(my.mod)$coefficients
      est.vector[i] = gee.values['trt:period', 'Estimate']
      se.vector[i] = gee.values['trt:period', 'Std.err']
      stat.vector[i] = gee.values['trt:period', 'Wald']
      pval.vector[i] = gee.values['trt:period', 'Pr(>|W|)']
      converge[i] <- ifelse(summary(my.mod)$error == 0, TRUE, FALSE)
    }
    
    # option to stop the function early if fits are singular
    if (poorFitOverride == FALSE & converge[i] == FALSE & i > 50) {
      if (sum(converge == FALSE, na.rm = TRUE) > (nsim * .25)) {
        stop("more than 25% of simulations are singular fit: check model specifications")
      }
    }
    
    # stop the loop if power is <0.5
    if (lowPowerOverride == FALSE & i > 50 & (i %% 10 == 0)) {
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
    " Simulations: Difference in Difference, Continuous Outcome"
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
    converge = converge
  )
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  cps.model.temp <- dplyr::filter(cps.model.est, converge == TRUE)
  power.parms <- confintCalc(alpha = alpha,
                             p.val = cps.model.temp[, 'p.value'])
  
  # Create object containing arm & time-specific differences
  values.vector = values.vector / nsim
  differences = data.frame(
    Period = c(0, 0, 1, 1),
    Arm.2 = c(0, 1, 0, 1),
    Values = round(values.vector, 3)
  )
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Arm.1' = nsubjects[1:nclusters[1]],
                       'Arm.2' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm.1" = c("n.clust" = nclusters[1]),
    "Arm.2" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing variance parameters for each group at each time point
  var.parms = list(
    "Time.Point.0" = data.frame(
      'Arm.1' = c("sigma_sq" = sigma_sq[1], "sigma_b_sq" = sigma_b_sq0[1]),
      'Arm.2' = c("sigma_sq" = sigma_sq[2], "sigma_b_sq" = sigma_b_sq0[2])
    ),
    "Time.Point.1" = data.frame(
      'Arm.1' = c("sigma_sq" = sigma_sq[3], "sigma_b_sq" = sigma_b_sq1[1]),
      'Arm.2' = c("sigma_sq" = sigma_sq[4], "sigma_b_sq" = sigma_b_sq1[2])
    )
  )
  
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
      "means" = c(mu, mu2),
      "model.estimates" = cps.model.est,
      "sim.data" = simulated.datasets,
      "differences" = differences,
      "convergence" = converge
    ),
    class = 'crtpwr'
  )
  return(complete.output)
}
