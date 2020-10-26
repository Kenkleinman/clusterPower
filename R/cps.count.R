#' Simulation-based power estimation for cluster-randomized trials: Parallel Designs, Count Outcome
#'
#' @description 
#' \loadmathjax
#' 
#'  
#' This function uses Monte Carlo methods (simulations) to estimate 
#' power for cluster-randomized trials with integer-valued outcomes. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, between-cluster variance, and
#' two of the following three parameters: mean event rate per unit time in one group,
#' the mean event rate per unit time in the second group, and/or the
#' mean difference in event rates between groups. Default values are provided
#' for significance level, analytic method, whether progress updates are displayed, 
#' and whether the simulated data sets are retained.
#'
#' Note that if all units have the same observation time, you can use the
#' mean count instead of the "mean event per unit time" in the preceding paragraph.
#'
#'
#'
#' @param nsim Number of datasets to simulate; accepts integer. Required.
#' 
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (implying equal cluster sizes for the two groups), 
#' a vector of length two (equal cluster sizes within arm), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within arm).  If a vector of > 2 is provided in
#' \code{nsubjects}, \code{sum(nclusters)} must match the \code{nsubjects} vector length.  Required.
#' 
#' @param nclusters Number of clusters per treatment group; accepts a single integer
#' (if there are the same number of clusters in each arm) or a vector of 2 integers
#' (if there are not). 
#' Required.
#' 
#' At least 2 of the following 3 arguments must be specified:
#' 
#' @param c1 The mean event rate per unit time in the first arm.
#' 
#' @param c2 The mean event rate per unit time in the second arm.
#' 
#' @param cDiff Expected difference in mean event rates between groups, defined as
#' \code{cDiff = c1 - c2}.
#' 
#' @param sigma_b_sq Between-cluster variance; if sigma_b_sq2 is not specified,
#' between-cluster variances are assumed to be equal in the two arms. Accepts numeric. Required.
#' 
#' @param sigma_b_sq2 Between-cluster variance for clusters in the second arm. Only required if 
#' between-cluster variances differ between treatment arms.
#' 
#' @param family Distribution from which responses are simulated. Accepts Poisson
#' ('poisson') or negative binomial ('neg.binom'); default = 'poisson'. Required.
#' 
#' @param negBinomSize Only used when generating simulated data from the 
#' negative binomial (family = 'neg.binom'), this is the target for number of 
#' successful trials, or the dispersion parameter (the shape parameter of the gamma 
#' mixing distribution). Must be strictly positive but need not be integer. 
#' Defaults to 1.
#' 
#' @param method Data analysis method, either generalized linear mixed effects model
#' (GLMM) or generalized estimating equations (GEE). Accepts c('glmm', 'gee'); default = 'glmm'.
#' Required.
#' 
#' @param analysis Family used for data analysis; currently only applicable when \code{method = 'glmm'}.
#' Accepts c('poisson', 'neg.binom'); default = 'poisson'. Required.
#' 
#' @param alpha The level of significance of the test, the probability of a
#' Type I error. Default = 0.05.
#' 
#' @param quiet When set to \code{FALSE}, displays simulation progress and estimated
#' completion time. Default = \code{FALSE}.
#' 
#' 
#' @param allSimData Option to include a list of all simulated datasets in the output object.
#' Default = \code{FALSE}.
#' 
#' @param seed Option to set the seed. Default is NA.
#' 
#' @param nofit Option to skip model fitting and analysis and instead return a dataframe with
#' the simulated datasets. Default = \code{FALSE}.
#' 
#' @param poorFitOverride Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge.
#' 
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' 
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' 
#' @param optimizer Option to fit with a different optimizer from the package
#' \code{optimx}. Defaults to L-BFGS-B. See optimx package documentation for all options.
#' 
#' @param irgtt Logical. Default = FALSE. Is the experimental design an 
#' individually randomized group treatment trial? For details, 
#' see ?cps.irgtt.count.
#' 
#'
#'
#'
#'
#' @return If \code{nofit = F}, a list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, distribution of
#'   simulated data, and regression family
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power),
#'                "lower.95.ci" (Lower 95% confidence interval bound),
#'                "upper.95.ci" (Upper 95% confidence interval bound).
#'                Note that non-convergent models are returned for review,
#'                but not included in this calculation.
#'   \item Analytic method used for power estimation
#'   \item Data frame containing families for distribution and analysis of simulated data
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting \mjseqn{\sigma_b^2} for each group
#'   \item Vector containing expected events per unit time and risk ratios based on user inputs
#'   \item Data frame with columns:
#'                   "Estimate" (Estimate of treatment effect for a given simulation),
#'                   "Std.err" (Standard error for treatment effect estimate),
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)),
#'                   "p.value",
#'                   "converge" (Did model converge for that set of simulated data?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing:
#'                   "y" (Simulated response value),
#'                   "trt" (Indicator for treatment arm),
#'                   "clust" (Indicator for cluster)
#'   \item Logical vector reporting whether models converged.
#' }
#' 
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' 
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#'   
#' @details 
#' 
#' 
#' If \code{family = 'poisson'}, the data generating model is:
#' \mjsdeqn{y_{ij} \sim \text{Poisson}(e^{c_1 + b_i}) }
#' for the first group or arm, where \mjseqn{b_i \sim N(0,\sigma_b^2)}, 
#' while for the second group, 
#'  
#' \mjsdeqn{y_{ij} \sim \text{Poisson}(e^{c_2 + b_i}) }
#' where \mjseqn{b_i \sim N(0,\sigma_{b_2}^2)}; if 
#' \mjseqn{\sigma_{b_2}^2} is not specified, then the second group uses
#' \mjseqn{b_i \sim N(0,\sigma_b^2)}.
#' 
#' If \code{family = 'neg.bin'}, the data generating model, using the
#' alternative parameterization of the negative binomial distribution
#' detailed in \code{stats::rnbinom}, is:
#' 
#' \mjsdeqn{y_{ij} \sim \text{NB}(\mu = e^{c_1 + b_i}, \text{size} = 1) }
#' 
#' for the first group or arm, where \mjseqn{b_i \sim N(0,\sigma_b^2)}, 
#' while for the second group, 
#'  
#' \mjsdeqn{y_{ij} \sim \text{NB}(\mu = e^{c_2 + b_i}, \text{size} = 1) }
#' where \mjseqn{b_i \sim N(0,\sigma_{b_2}^2)}; if 
#' \mjseqn{\sigma_{b_2}^2} is not specified, then the second group uses
#' \mjseqn{b_i \sim N(0,\sigma_b^2)}.
#' 
#' 
#' 
#' 
#' Non-convergent models are not included in the calculation of exact confidence 
#' intervals.
#' 
#' 
#' 
#' @section Testing details:
#' This function has been verified, where possible, against reference values from PASS11,
#' \code{CRTsize::n4incidence}, \code{clusterPower::cps.count}, and
#' \code{clusterPower::cpa.count}.
#' 
#' @author Alexander R. Bogdan, Alexandria C. Sakrejda 
#' (\email{acbro0@@umass.edu}), and Ken Kleinman 
#' (\email{ken.kleinman@@gmail.com})
#'
#'
#' @examples
#' 
#' # Estimate power for a trial with 10 clusters in each arm with 20 subjects each, 
#' # with sigma_b_sq = 0.1 in both arms. We expect mean event rates per unit time of
#' # 20 and 30 in the first and second arms, respectively, and we use 100 simulated
#' # data sets analyzed by the GEE method. 
#' 
#' \dontrun{
#' count.sim = cps.count(nsim = 100, nsubjects = 20, nclusters = 10,
#'                       c1 = 20, c2 = 30, sigma_b_sq = 0.1,
#'                       family = 'poisson', analysis = 'poisson',
#'                       method = 'gee', alpha = 0.05, quiet = FALSE,
#'                       allSimData = FALSE, seed = 123)
#' }                  
#' # The resulting estimated power (if you set seed = 123) should be about 0.8.
#'
#'
#'
#' # Estimate power for a trial with 10 clusters and 10 subjects per cluster in the
#' # first arm, 20 clusters and 20 subjects per cluster in the second, and
#' # sigma_b_sq = 0.1 in both arms. We expect mean event rates per unit time of
#' # 20 and 30 in the first and second arms, respectively, and we use 100 simulated
#' # data sets analyzed by the GLMM method. 
#' 
#' \dontrun{
#' count.sim = cps.count(nsim = 100, nsubjects = c(10,20), nclusters = c(10,20),
#'                       c1 = 20, c2 = 30, sigma_b_sq = 0.1,
#'                       family = 'poisson', analysis = 'poisson',
#'                       method = 'glmm', alpha = 0.05, quiet = FALSE,
#'                       allSimData = FALSE, seed = 123)
#' }               
#' # The resulting estimated power (if you set seed = 123) should be about 0.95.
#'
#'  
#' 
#' # Estimate power for a trial with 5 clusters in the first arm, those clusters having
#' # 4, 5, 6, 7, and 7 subjects each, and 10 clusters in the second arm, those
#' # clusters having 5 subjects each, with sigma_b_sq = 0.1 in the first arm and
#' # sigma_b_sq2 = .05 in the second arm. We expect mean event rates per unit time
#' # of 20 and 30 in the first and second arms, respectively, and we use 100 simulated
#' # data sets analyzed by the GLMM method.                      
#' 
#' \dontrun{                                            
#' count.sim = cps.count(nsim = 100, nsubjects = c(4, 5, 6, 7, 7, rep(5, times = 10)), 
#'                       nclusters = c(5,10),
#'                       c1 = 20, c2 = 30,
#'                       sigma_b_sq = 0.1, sigma_b_sq2 = 0.05,
#'                       family = 'poisson', analysis = 'poisson',
#'                       method = 'glmm', alpha = 0.05, quiet = FALSE,
#'                       allSimData = FALSE, seed = 123, optimizer = "L-BFGS-B")
#' }
#' # The resulting estimated power (if you set seed = 123) should be about 0.75.
#'
#'
#'
#' @export

# Define function
cps.count = function(nsim = NULL,
                     nsubjects = NULL,
                     nclusters = NULL,
                     c1 = NULL,
                     c2 = NULL,
                     cDiff = NULL,
                     sigma_b_sq = NULL,
                     sigma_b_sq2 = NULL,
                     family = 'poisson',
                     negBinomSize = 1,
                     analysis = 'poisson',
                     method = 'glmm',
                     alpha = 0.05,
                     quiet = FALSE,
                     allSimData = FALSE,
                     irgtt = FALSE,
                     seed = NA,
                     nofit = FALSE,
                     poorFitOverride = FALSE,
                     lowPowerOverride = FALSE,
                     timelimitOverride = TRUE,
                     optimizer = "L-BFGS-B") {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  # Create vectors to collect iteration-specific values
  est.vector <- NULL
  se.vector <- NULL
  stat.vector <- NULL
  pval.vector <- NULL
  converge.vector <- NULL
  simulated.datasets <- list()
  converge <- NULL
  
  # Create progress bar
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
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS
  sim.data.arg.list = list(nsim, nclusters, nsubjects, sigma_b_sq)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if (sum(sim.data.args) > 0) {
    stop(
      "NSIM, NSUBJECTS, NCLUSTERS & sigma_b_sq must all be specified. Please review your input values."
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
  # Set cluster sizes for arm (if not already specified)
  if (length(nclusters) == 1) {
    if (irgtt == TRUE) {
      nclusters[2] = nclusters[1]
      nclusters[1] = 1
    } else {
      nclusters[2] = nclusters[1]
    }
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
  
  # Validate sigma_b_sq, sigma_b_sq2, ALPHA
  if (length(sigma_b_sq) > 1 || length(sigma_b_sq2) > 1) {
    errorCondition("The lengths of sigma_b_sq and sigma_b_sq2 cannot be larger than 1.")
  }
  if (irgtt == FALSE) {
    min0.warning = " must be a numeric value greater than 0"
    if (!is.numeric(sigma_b_sq) || sigma_b_sq <= 0) {
      stop("sigma_b_sq", min0.warning)
    }
    if (!is.null(sigma_b_sq2) && sigma_b_sq2 <= 0) {
      stop("sigma_b_sq2", min0.warning)
    }
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate C1, C2, cDiff
  parm1.arg.list = list(c1, c2, cDiff)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if (sum(parm1.args) > 1) {
    stop("At least two of the following terms must be specified: C1, C2, cDiff")
  }
  if (sum(parm1.args) == 0 && cDiff != abs(c1 - c2)) {
    stop("At least one of the following terms has be misspecified: C1, C2, cDiff")
  }
  
  # Validate FAMILY, ANALYSIS, METHOD, QUIET
  if (!is.element(family, c('poisson', 'neg.binom'))) {
    stop(
      "FAMILY must be either 'poisson' (Poisson distribution)
         or 'neg.binom'(Negative binomial distribution)"
    )
  }
  if (!is.element(analysis, c('poisson', 'neg.binom'))) {
    stop(
      "ANALYSIS must be either 'poisson' (Poisson regression)
         or 'neg.binom'(Negative binomial regression)"
    )
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
  if (family == 'neg.binom' && negBinomSize < 0) {
    stop(
      "negBinomSize must be positive"
    )
  }
  
  # Calculate inputs & variance parameters
  if (is.null(c1)) {
    c1 = abs(cDiff - c2)
  }
  if (is.null(c2)) {
    c2 = abs(c1 - cDiff)
  }
  if (is.null(cDiff)) {
    cDiff = c1 - c2
  }
  if (is.null(sigma_b_sq2)) {
    sigma_b_sq[2] = sigma_b_sq
  } else{
    sigma_b_sq[2] = sigma_b_sq2
  }
  
  # Create indicators for arm & cluster
  trt = c(rep(1, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(2, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  
  est.vector <- vector(mode = "numeric", length = nsim)
  se.vector <- vector(mode = "numeric", length = nsim)
  stat.vector <- vector(mode = "numeric", length = nsim)
  pval.vector <- vector(mode = "numeric", length = nsim)
  converge.vector <- vector(mode = "logical", length = nsim)

  # Create simulation loop
  for (i in 1:nsim) {
    # Generate between-cluster effects for arm 1 and arm 2
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq[2]))

    # Create arm 1 y-value
    y0.intercept <- list()
    for (j in 1:length(nsubjects[1:nclusters[1]])) {
      y0.intercept[[j]] <- rep(randint.0[j], each = nsubjects[j])
    }
    y0.intercept <- unlist(y0.intercept)
    y0.linpred = y0.intercept + log(c1)
    y0.prob = exp(y0.linpred)
    if (family == 'poisson') {
      y0 = stats::rpois(length(y0.prob), y0.prob)
    }
    if (family == 'neg.binom') {
      y0 = stats::rnbinom(length(y0.prob), size = negBinomSize, mu = y0.prob)
    }

    # Create arm 2 y-value
    y1.intercept <- list()
    for (j in 1:length(nsubjects[(nclusters[1] + 1):length(nsubjects)])) {
      y1.intercept[[j]] <- rep(randint.1[j], each = nsubjects[nclusters[1] + j])
    }
    y1.intercept <- unlist(y1.intercept)
    y1.linpred = y1.intercept + log(c2) #+ log((c1 / (1 - c1)) / (c2 / (1 - c2)))
    y1.prob = exp(y1.linpred)
    if (family == 'poisson') {
      y1 = stats::rpois(length(y1.prob), y1.prob)
    }
    if (family == 'neg.binom') {
      y1 = stats::rnbinom(length(y1.prob), size = 1, mu = y1.prob)
    }
    
    # Create single response vector
    y = c(y0, y1)

    # Create and store data for simulated dataset
    sim.dat = data.frame(trt = as.factor(trt),
                         clust = as.factor(clust),
                         y = as.integer(y))
    if (allSimData == TRUE) {
      simulated.datasets[[i]] = list(sim.dat)
    }
    
    # option to return simulated data only
    if (nofit == TRUE) {
      if (!exists("nofitop")) {
        nofitop <- data.frame(trt = trt,
                              clust = clust,
                              y1 = y)
      } else {
        nofitop[, length(nofitop) + 1] <- y
      }
      if (length(nofitop) == (nsim + 2)) {
        temp1 <- seq(1:nsim)
        temp2 <- paste0("y", temp1)
        colnames(nofitop) <- c("arm", "cluster", temp2)
      }
      if (length(nofitop) != (nsim + 2)) {
        next()
      }
      return(nofitop)
    }
    
    # Set warnings to OFF
    # Note: Warnings will still be stored in 'warning.list'
    options(warn = -1)
    
    #set start time
    start.time = Sys.time()
    
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      require("lme4")
      if (i == 1) {
        require("optimx")
        if (isTRUE(optimizer == "auto")) {
          if (irgtt == FALSE) {
            if (analysis == 'poisson') {
              my.mod = try(lme4::glmer(
                y ~ as.factor(trt) + (1 | clust),
                data = sim.dat,
                family = stats::poisson(link = 'log')
              ))
            }
            if (analysis == 'neg.binom') {
              my.mod = try(lme4::glmer.nb(y ~ as.factor(trt) + (1 |
                                                              clust), data = sim.dat))
            }
          }
          if (irgtt == TRUE) {
            if (analysis == 'poisson') {
              my.mod <-
                try(lme4::glmer(
                  y ~ trt + (0 + as.factor(trt) | clust),
                  data = sim.dat,
                  family = stats::poisson(link = 'log'))
                )
            }
            if (analysis == 'neg.binom') {
              my.mod = try(lme4::glmer.nb(y ~ trt + (0 + as.factor(trt) |
                                                   clust), data = sim.dat))
            }
          }
          goodopt <- optimizerSearch(my.mod)
        } else {
          goodopt <- optimizer
        }
      }
      if (irgtt == FALSE) {
        if (analysis == 'poisson') {
          my.mod = try(lme4::glmer(
            y ~ trt + (1 | clust),
            data = sim.dat,
            family = stats::poisson(link = 'log'),
            control = lme4::glmerControl(
              optimizer = "optimx",
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
          )
        }

        if (analysis == 'neg.binom') {
          my.mod = try(lme4::glmer.nb(
            y ~ trt + (1 | clust),
            data = sim.dat,
            control = lme4::glmerControl(
              optimizer = "optimx",
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
          )
        }
      }
      if (irgtt == TRUE) {
        if (analysis == 'poisson') {
          my.mod <- try(lme4::glmer(
            y ~ trt + (0 + trt | clust),
            data = sim.dat,
            family = stats::poisson(link = 'log'),
            control = lme4::glmerControl(
              optimizer = "optimx",
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
          )
        }
        if (analysis == 'neg.binom') {
          my.mod = try(lme4::glmer.nb(
            y ~ trt + (0 + trt | clust),
            data = sim.dat,
            control = lme4::glmerControl(
              optimizer = "optimx",
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
          )
        }
      }
      if (class(my.mod) == "try-error") {
        next
      }
      glmm.values <- summary(my.mod)$coefficient
      est.vector[i] <- glmm.values['trt2', 'Estimate']
      se.vector[i] <- glmm.values['trt2', 'Std. Error']
      stat.vector[i] <- glmm.values['trt2', 'z value']
      pval.vector[i] <- glmm.values['trt2', 'Pr(>|z|)']
      converge.vector[i] <- ifelse(
        is.null(my.mod@optinfo$conv$lme4$messages), TRUE, FALSE)
    }
    # Fit GEE (geeglm)
    if (method == 'gee') {
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(
        y ~ trt,
        data = sim.dat,
        family = stats::poisson(link = 'log'),
        id = clust,
        corstr = "exchangeable"
      )
      gee.values <- summary(my.mod)$coefficients
      est.vector[i] <- gee.values['trt2', 'Estimate']
      se.vector[i] <- gee.values['trt2', 'Std.err']
      stat.vector[i] <- gee.values['trt2', 'Wald']
      pval.vector[i] <- gee.values['trt2', 'Pr(>|W|)']
      converge.vector[i] <- ifelse(my.mod$geese$error != 0, FALSE, TRUE)
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
    
    # option to stop the function early if fits are singular
    if (poorFitOverride == FALSE && converge.vector[i] == FALSE && i > 12) {
      if (sum(converge.vector == FALSE, na.rm = TRUE) > (nsim * .25)) {
        stop(
          "more than 25% of simulations are singular fit: check model specifications"
        )
      }
    }

    # Set warnings to ON
    # Note: Warnings stored in 'warning.list'
    options(warn = 0)
    
    # Update simulation progress information
    if (quiet == FALSE) {
      if (i == 1) {
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
        
        #time limit override (for Shiny)
        if (min.est > 2 && timelimitOverride == FALSE){
          stop(paste0("Estimated completion time: ",
                      hr.est,
                      'Hr:',
                      min.est,
                      'Min'
          ))
        }
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
      # Iterate progress bar
      prog.bar$update(i / nsim)
      Sys.sleep(1 / 100)
      
      if (i == nsim) {
        message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
      }
    }
  }
  
  ## Output objects
  # Create object containing summary statement
  if (irgtt == FALSE) {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: Simple Design, Count Outcome. Data Simulated from ",
      switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
      " distribution. Analyzed using ",
      switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
      " regression"
    )
  } else {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: IRGTT Design, Count Outcome. Data Simulated from ",
      switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
      " distribution. Analyzed using ",
      switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
      " regression"
    )
  }
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model',
                       gee = 'Generalized Estimating Equation')
  
  # Store simulation output in data frame
  cps.model.est = data.frame(
    Estimate = as.vector(unlist(est.vector)),
    Std.err = as.vector(unlist(se.vector)),
    Test.statistic = as.vector(unlist(stat.vector)),
    p.value = as.vector(unlist(pval.vector)),
    converge = as.vector(unlist(converge.vector))
  )
  
  # Calculate and store power estimate & confidence intervals
  cps.model.temp <- dplyr::filter(cps.model.est, converge == TRUE)

  power.parms <- confintCalc(alpha = alpha,
                              p.val = cps.model.temp[, 'p.value'])
  
  # Create object containing inputs
  c1.c2.rr = round(exp(log(c1) - log(c2)), 3)
  c2.c1.rr = round(exp(log(c2) - log(c1)), 3)
  inputs = t(data.frame(
    'Arm1' = c("count" = c1, "risk.ratio" = c1.c2.rr),
    'Arm2' = c("count" = c2, 'risk.ratio' = c2.c1.rr),
    'Difference' = c("count" = cDiff, 'risk.ratio' = c2.c1.rr - c1.c2.rr)
  ))

  # Create object containing group-specific cluster sizes
  cluster.sizes = list(
    'Arm1' = nsubjects[1:nclusters[1]],
    'Arm2' = nsubjects[(1 + nclusters[1]):(nclusters[1] + nclusters[2])]
  )
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm1" = c("n.clust" = nclusters[1]),
    "Arm2" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame(
    'Arm1' = c('sigma_b_sq' = sigma_b_sq[1]),
    'Arm2' = c('sigma_b_sq' = sigma_b_sq[2])
  ))
  
  # Create object containing FAMILY & REGRESSION parameters
  dist.parms = rbind(
    'Family:' = paste0(switch(
      family, poisson = 'Poisson', neg.binom = 'Negative Binomial'
    ), ' distribution'),
    'Analysis:' = paste0(switch(
      analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'
    ), ' distribution')
  )
  colnames(dist.parms) = "Data Simuation & Analysis Parameters"
  
  
  # Create list containing all output (class 'crtpwr') and return
  complete.output = structure(
    list(
      "overview" = summary.message,
      "nsim" = nsim,
      "power" = power.parms,
      "method" = long.method,
      "dist.parms" = dist.parms,
      "alpha" = alpha,
      "cluster.sizes" = cluster.sizes,
      "n.clusters" = n.clusters,
      "variance.parms" = var.parms,
      "inputs" = inputs,
      "model.estimates" = cps.model.est,
      "sim.data" = simulated.datasets,
      "convergence" = converge.vector
    ),
    class = 'crtpwr'
  )
  return(complete.output)
}
