#' Power simulations for cluster-randomized trials: Parallel Designs, Normal Outcome
#' 
#' @description 
#' \loadmathjax
#' 
#' This function uses Monte Carlo methods (simulations) to estimate 
#' power for parallel design cluster-randomized trials with normal outcomes. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per arm, expected means of the arms, and two of 
#' the following: ICC, within-cluster variance, or between-cluster variance.  
#' Defaults are provided for significance level, analytic method, progress 
#' updates, and whether the simulated data sets are retained.
#' 
#' Users have the option of specifying different variance parameters for each
#' arm, different numbers of clusters for each treatment group, and different numbers
#' of units within each cluster. 
#' 
#' Non-convergent models are not included in the calculation of exact confidence 
#' intervals.
#' 
#' @section Testing details:   
#' This function has been verified, where possible, against reference values from the NIH's GRT 
#' Sample Size Calculator, PASS11, \code{CRTsize::n4means}, and 
#' \code{clusterPower::cpa.normal}.
#' 
#' @param nsim Number of datasets to simulate; accepts integer.  Required.
#' 
#' @param nclusters Number of clusters per condition; accepts single integer (implying equal numbers of clusters in the two groups)
#' or vector of length 2 (unequal number of clusters per arm).  Required.
#' 
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (implying equal cluster sizes for the two groups), 
#' a vector of length two (equal cluster sizes within arm), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within arm).  Required.
#' 
#' @param mu Mean in the first arm; accepts numeric, default 0.  Required..
#' 
#' @param mu2 Mean in the second arm; accepts numeric.  Required.
#' 
#' At least 2 of the following must be specified:
#' 
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 and 1.
#' 
#' @param sigma_sq Within-cluster variance; accepts numeric.
#' 
#' @param sigma_b_sq Between-cluster variance; accepts numeric.
#' 
#' 
#' The defaults for the following are all NA, implying equal variance parameters 
#' for the two groups. If one of the following is given, variance parameters differ 
#' between treatment groups, and at least 2 of the following 
#' must be specified:
#' 
#' @param ICC2 Intra-cluster correlation coefficient for clusters in the second arm.
#' 
#' @param sigma_sq2 Within-cluster variance for clusters in the second arm.
#' 
#' @param sigma_b_sq2 Between-cluster variance for clusters in the second arm.
#' 
#' Optional parameters:
#' 
#' @param alpha Significance level; default = 0.05.
#' 
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM, default) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee').
#' 
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' 
#' @param allSimData Option to include a list of all simulated datasets in the output object.
#' Default = \code{FALSE}.
#' 
#' @param seed Option to set the seed. Default, NA, selects a seed based on the system clock.
#' 
#' @param irgtt Logical. Is the experimental design an individually randomized 
#' group treatment trial? For details, see ?cps.irgtt.normal.
#' 
#' @param poorFitOverride Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge.
#' 
#' @param nofit Option to skip model fitting and analysis and instead return a dataframe with
#' the simulated datasets. Default = \code{FALSE}.
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
#' @return If \code{nofit = F}, a list with the following components:
#' \itemize{
#'   \item Character string indicating total number of simulations and simulation type
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95\% confidence interval bound), 
#'                "upper.95.ci" (Upper 95\% confidence interval bound),
#'                "Alpha" (Probability of committing a type I or \mjseqn{\alpha} error or rejecting a true null),
#'                "Beta" (Probability of committing a type II error or failing to reject a false null).
#'                Note that non-convergent models are returned for review, 
#'                but not included in this calculation.
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters in each arm
#'   \item Data frame reporting ICC, variance parameters, and means for each arm
#'   \item Vector containing expected group means based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "converge", (Did the model converge?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for arm), 
#'                   "clust" (Indicator for cluster)
#'                   }
#' 
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "clust" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#' 
#' @details 
#'
#' The data generating model for observation \mjseqn{i} in cluster \mjseqn{j} is:
#' \mjsdeqn{y_{ij} \sim N(\mu + b_i, \sigma^2) }
#' for the first group or arm, where \mjseqn{b_i \sim N(0,\sigma_b^2)} 
#' , while for the second group, 
#'  
#' \mjsdeqn{y_{ij} \sim N(\mu_2 + b_i, \sigma_2^2) }
#' where \mjseqn{b_i \sim N(0,\sigma_{b_2}^2)}; if none of 
#' \mjseqn{\sigma_2^2, \sigma_{b_2}^2} or \code{ICC2} are used, then the second group uses
#' \mjseqn{b_i \sim N(0,\sigma_b^2)} 
#' and \mjseqn{y_{ij} \sim N(\mu_2 + b_i, \sigma^2)} .
#' 
#' All random terms are generated indepedent of one another.
#' 
#' 
#' For calls without \mjseqn{\sigma_2^2, \sigma_{b_2}^2} or \code{ICC2}, and using
#' \code{method="glmm"} the fitted model is:
#' \mjsdeqn{y_{ij}|b_i = \mu + \beta_1 x_{ij}  + b_i + e_{ij}}
#'
#' with \mjseqn{\beta_1 = \mu_2 - \mu},
#' treatment group indicator \mjseqn{x_{ij} = 0} for the first group, with
#' \mjseqn{b_i \sim N(0, \sigma_b^2)} and \mjseqn{e_{ij} \sim N(0,\sigma^2)}.
#' In this case, both the random effects distribution and the residual distribution are the same for both
#' conditions.
#' 
#' Otherwise, for \code{method="glmm"} the fitted model is:
#' \mjsdeqn{y_{ij}|b_i = \mu + \beta_1 x_{ij} 
#'    + b_i I(x_{ij}=0) + e_{ij} I(x_{ij}=0) 
#'    + g_i I(x_{ij}=1) + f_{ij} I(x_{ij}=1)
#'    }
#'
#' with \mjseqn{\beta_1}, \mjseqn{x_{ij}, b_i}, and \mjseqn{e_{ij}} as above, with
#' \mjseqn{g_i \sim N(0, \sigma_{b_2}^2)} and \mjseqn{f \sim N(0,\sigma_2^2)}, the 
#' random effects and residual distribution in the second group.
#' 
#' @examples
#' 
#' # Estimate power for a trial with 10 clusters in each arm and 25 subjects in each 
#' # cluster, with an ICC of .3, sigma squared of 20 (implying sigma_b^2 of 8.57143) 
#' # in each group, with arm means of 1 and 4.75 in the two groups, using 100 simulated 
#' # data sets. The resulting estimated power should be 0.78.
#'    
#' \dontrun{
#' 
#' normal.sim = cps.normal(nsim = 100, nsubjects = 25, nclusters = 10, mu = 1, 
#'   mu2 = 4.75, ICC = 0.3, sigma_sq = 20, seed = 123)
#' }
#' 
#' 
#' 
#' # Estimate power for a trial with 5 clusters in one arm, those clusters having 25 subjects 
#' # each, 25 clusters in the other arm, those clusters having 5 subjects each, the first arm
#' # having a sigma squared of 20 and sigma_b squared of 8.57143, and the second a sigma squared
#' # of 9 and a sigma_b squared of 1, with estimated arm means of 1 and 4.75 in the first and 
#' # second groups, respectively, using 100 simulated data sets analyzed by the GEE method. 
#' # The estimated power should be 0.79, assuming seed = 123.
#' 
#' \dontrun{
#' normal.sim2 = cps.normal(nsim = 100, nclusters = c(5,25), nsubjects = c(25,5), mu = 1, 
#'   mu2 = 4.75, sigma_sq = 20,sigma_b_sq = 8.8571429, sigma_sq2 = 9, sigma_b_sq2 = 1, 
#'   method = "gee", seed = 123)
#' }
#' 
#' 
#' # Estimate power for a trial with 5 clusters in one arm, those clusters having
#' # 4, 5, 6, 7, 7, and 7 subjects each, and 10 clusters in the other arm,
#' # those clusters having 5 subjects each, with sigma_b_sq = .3 and and ICC of .3 in both arms.
#' # We have estimated arm means of 1 and 2 in the first and second arms, respectively, and we use
#' # 100 simulated data sets analyzed by the GLMM method.
#' 
#' \dontrun{
#' normal.sim2 = cps.normal(nsim = 100, nclusters = c(6,10), 
#'   nsubjects = list(c(4, 5, 6, 7, 7, 7), rep(5, times = 10)),
#'   mu = 1, mu2 = 2, sigma_b_sq = .3, ICC = .3, method = "glmm",
#'   seed = 1)
#' }
#' 
#' # The resulting estimated power (if you set seed = 1) should be about 0.76.
#' 
#' # Estimate power for a trial with 3 clusters in one arm, those clusters having 25, 35, and 45 subjects 
#' # each, and 10 clusters in the other arm, those clusters having 5 subjects each, the first arm
#' # having a sigma squared of 20 and sigma_b squared of 8.57143, and the second a sigma squared
#' # of 9 and a sigma_b squared of 1, with estimated arm means of 1 and 4.75 in the first and 
#' # second groups, respectively, using 100 simulated data sets analyzed by the GLMM method.
#' 
#' \dontrun{
#' 
#' normal.sim2 = cps.normal(nsim = 100, nclusters = c(3,10), 
#'   nsubjects = c(25, 35, 45, rep(5, times = 10)),
#'   mu = 1, mu2 = 4.75, sigma_sq = 20, sigma_b_sq = 8.8571429,
#'   sigma_sq2 = 9, sigma_b_sq2 = 1, method = "glmm")
#' }
#' 
#' # The resulting estimated power (if you set seed = 1) should be about 0.71.
#' 
#' 
#' @author Alexander R. Bogdan, Alexandria C. Sakrejda 
#' (\email{acbro0@@umass.edu}), and Ken Kleinman 
#' (\email{ken.kleinman@@gmail.com})
#' 
#' 
#' 
#' 
#' 
#' @export
cps.normal = function(nsim = NA,
                      nclusters = NA,
                      nsubjects = NA,
                      mu = 0,
                      mu2 = NA,
                      ICC = NA,
                      sigma_sq = NA,
                      sigma_b_sq = NA,
                      ICC2 = NA,
                      sigma_sq2 = NA,
                      sigma_b_sq2 = NA,
                      alpha = 0.05,
                      method = 'glmm',
                      quiet = FALSE,
                      allSimData = FALSE,
                      seed = NA,
                      poorFitOverride = FALSE,
                      timelimitOverride = TRUE,
                      lowPowerOverride = FALSE,
                      irgtt = FALSE,
                      nofit = FALSE) {
  # option for reproducibility
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create vectors to collect iteration-specific values
  est.vector = rep(NA, length = nsim)
  se.vector = rep(NA, length = nsim)
  stat.vector = rep(NA, length = nsim)
  pval.vector = rep(NA, length = nsim)
  # This container keeps track of how many models failed to converge
  converge.vector <- rep(NA, length = nsim)
  simulated.datasets = list()
  
  # initialize progress bar
  prog.bar =  progress::progress_bar$new(
    format = "(:spin) [:bar] :percent eta :eta",
    total = nsim,
    clear = FALSE,
    width = 100
  )
  prog.bar$tick(0)
  
  # Validate NSIM, NCLUSTERS, NSUBJECTS
  sim.data.arg.list = list(nsim, nclusters, nsubjects)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.na))
  if (sum(sim.data.args) > 0) {
    stop("NSIM, NCLUSTERS, & NSUBJECTS must all be specified. Please review your input values.")
  }
  min1.warning = " must be an integer greater than or equal to 1"
  if (!is.wholenumber(nsim) || nsim < 1) {
    stop(paste0("NSIM", min1.warning))
  }

  if (!is.wholenumber(nclusters) || nclusters < 1) {
    stop(paste0("NCLUSTERS", min1.warning))
  }
  if (is.list(nsubjects)){
    temp <- unlist(nsubjects)
  } else {
    temp <- nsubjects
  }
  if (!is.wholenumber(temp) || temp < 1) {
    stop(paste0("NSUBJECTS", min1.warning))
  }
  if (length(nclusters) > 2) {
    stop(
      "NCLUSTERS can only be a scalar (equal # of clusters per group) or a vector of length 2 (unequal # of clusters per group)"
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
  
  ## Create variance parameters
  # sigma_b_sq, sigma_sq, ICC
  if (!is.na(c(ICC, sigma_sq)) && is.na(sigma_b_sq)) {
    sigma_b_sq = ICC * sigma_sq / (1 - ICC)
  }
  if (!is.na(c(ICC, sigma_b_sq)) && is.na(sigma_sq)) {
    sigma_sq = sigma_b_sq / ICC - sigma_b_sq
  }

  if (!is.na(c(sigma_sq, sigma_b_sq)) && is.na(ICC)) {
    ICC = sigma_b_sq / (sigma_b_sq + sigma_sq)
  }
  # sigma_b_sq2, sigma_sq2, ICC2
  if (!is.na(c(ICC2, sigma_sq2)) && is.na(sigma_b_sq2)) {
    sigma_b_sq2 = ICC2 * sigma_sq2 / (1 - ICC2)
  }
  if (!is.na(c(ICC2, sigma_b_sq2)) && is.na(sigma_sq2)) {
    sigma_sq2 = sigma_b_sq2 / ICC2 - sigma_b_sq2
  }
  if (!is.na(c(sigma_sq2, sigma_b_sq2)) && is.na(ICC2)) {
    ICC2 = sigma_b_sq2 / (sigma_b_sq2 + sigma_sq2)
  }
  
  # Set within/between cluster variances & ICC for arm (if not already specified)
  if (isTRUE(is.na(sigma_sq2))) {
    sigma_sq2 <- sigma_sq
  }
  if (isTRUE(is.na(sigma_b_sq2))) {
    sigma_b_sq2 <- sigma_b_sq
  }
  if (isTRUE(is.na(ICC2))) {
    ICC2 <- ICC
  }
  
  # Validate mu, mu2, ALPHA
  if (is.na(mu) || is.na(mu2)) {
    stop("MU and MU2 are required.")
  }
  min0.warning = " must be numeric."
  if (!is.numeric(mu) || !is.numeric(mu2)) {
    stop("MU and MU2", min0.warning)
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate ICC, sigma_sq, sigma_b_sq, ICC2, sigma_sq2, sigma_b_sq2
  
  parm1.arg.list = list(ICC, sigma_sq, sigma_b_sq)
  parm1.args = unlist(lapply(parm1.arg.list, is.na))
  if (sum(parm1.args) > 1) {
    stop("At least two of the following terms must be specified: ICC, sigma_sq, sigma_b_sq")
  }

  if (round(ICC, 2) != round((sigma_b_sq / (sigma_b_sq + sigma_sq)), 2)) {
    stop("At least one of the following terms has been misspecified: ICC, sigma_sq, sigma_b_sq")
  }
  
  parm2.arg.list = list(ICC2, sigma_sq2, sigma_b_sq2)
  parm2.args = unlist(lapply(parm2.arg.list, is.na))
  if (sum(parm2.args) > 1 && sum(parm2.args) != 3) {
    stop(
      "At least two of the following terms must be provided to simulate arm-specific
         variances: ICC2, sigma_sq2, sigma_b_sq2"
    )
  }
  if (round(ICC2, 2) != round((sigma_b_sq2 / (sigma_b_sq2 + sigma_sq2)), 2)) {
    stop(
      "At least one of the following terms has been misspecified: ICC2, sigma_sq2, sigma_b_sq2"
    )
  }
  
  # Validate METHOD, QUIET, ALL.SIM.DATA
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
      "ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output"
    )
  }
  
  # Create indicators for arm & cluster
  if (is.list(nsubjects)){
    nsubjects <- unlist(nsubjects)
  }
  trt = c(rep(1, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(2, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))

  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  # Create simulation loop
  for (i in 1:nsim) {
    # Generate between-cluster effects
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq2))
    
    # Create y-value for the first arm
    y0.bclust = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.0[x], length.out = nsubjects[x])))
    y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x)
      stats::rnorm(
        x, mean = mu, sd = sqrt(sigma_sq)
      )))
    y.0 = y0.bclust + y0.wclust
    
    # Create y-value for the second arm
    y1.bclust = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                              function(x)
                                stats::rnorm(
                                  x, mean = mu2, sd = sqrt(sigma_sq2)
                                )))
    y.1 = y1.bclust + y1.wclust
    
    # Create single response vector
    y = c(y.0, y.1)

    # Create data frame for simulated dataset
    sim.dat = data.frame(y = y, trt = trt, clust = clust)
    if (allSimData == TRUE) {
      simulated.datasets[[i]] = sim.dat
    }
    
    # option to return simulated data only
    if (nofit == TRUE) {
      if (i == 1) {
        nofitop <- data.frame(trt = trt,
                              clust = clust,
                              y1 = y)
      } else {
        nofitop[, i + 2] <- y
      }
      if (i != nsim) {
        next()
      }
      if (i == nsim) {
        temp1 <- seq(1:nsim)
        temp2 <- paste0("y", temp1)
        colnames(nofitop) <- c("arm", "clust", temp2)
        return(nofitop)
      }
    }
    
    #set start time
    start.time = Sys.time()
    
    # trt and clust are re-coded as trt2 and clust2 to work nicely with lme.
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      if (irgtt == TRUE) {
        if (sigma_sq != sigma_sq2 && sigma_b_sq != sigma_b_sq2) {
          trt2 <- unlist(trt)
          clust2 <- unlist(clust)
          my.mod <-
            try(nlme::lme(
              y ~ as.factor(trt2),
              random =  ~ 0 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ),
            silent = TRUE)
          glmm.values <-  summary(my.mod)$tTable
          # get the overall p-values (>Chisq)
          null.mod <-
            try(nlme::lme(
              y ~ 1,
              random =  ~ 0 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ),
            silent = TRUE)

          pval.vector[i] = glmm.values['as.factor(trt2)2', 'p-value']
          est.vector[i] = glmm.values['as.factor(trt2)2', 'Value']
          se.vector[i] = glmm.values['as.factor(trt2)2', 'Std.Error']
          stat.vector[i] = glmm.values['as.factor(trt2)2', 't-value']
          converge.vector[i] <-
            ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE)
        }
        
        if (sigma_sq == sigma_sq2 && sigma_b_sq != sigma_b_sq2) {
          my.mod <-
            lmerTest::lmer(y ~ trt + (0 + as.factor(trt) | clust),
                           REML = FALSE,
                           data = sim.dat)
          # get the overall p-values (>Chisq)
          null.mod <-
            stats::update.formula(my.mod, y ~ 1 + (0 + as.factor(trt) |
                                                     clust))
          glmm.values[i] = summary(my.mod)$coefficients
          pval.vector[i] = glmm.values['trt', 'Pr(>|t|)']
          est.vector[i] = glmm.values['trt', 'Estimate']
          se.vector[i] = glmm.values['trt', 'Std. Error']
          stat.vector[i] = glmm.values['trt', 't value']
          converge.vector[i] = ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == FALSE, TRUE, FALSE)
          
          # option to stop the function early if fits are singular
          if (poorFitOverride == FALSE && converge.vector[i] == FALSE) {
            if (sum(converge.vector == FALSE, na.rm = TRUE) > (nsim * .25)) {
              stop(
                "more than 25% of simulations are singular fit: check model specifications"
              )
            }
            
          }
        }
        #if not IRGTT, then the following:
      } else {
        if (sigma_sq != sigma_sq2 && sigma_b_sq != sigma_b_sq2) {
          trt2 <- unlist(trt)
          clust2 <- unlist(clust)
          oldw <- getOption("warn")
          options(warn = -1)
          my.mod <-
            try(nlme::lme(
              y ~ as.factor(trt2),
              random =  ~ 1 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ),
            silent = TRUE)
          if (class(my.mod) != "try-error") {
            glmm.values <-  summary(my.mod)$tTable
            # get the overall p-values (>Chisq)
            null.mod <-
              try(nlme::lme(
                y ~ 1,
                random =  ~ 1 + as.factor(trt2) | clust2,
                weights = nlme::varIdent(form =  ~ 1 |
                                           as.factor(trt2)),
                method = "ML",
                control = nlme::lmeControl(opt = 'optim')
              ),
              silent = TRUE)
            options(warn = oldw)
            pval.vector[i] = glmm.values['as.factor(trt2)1', 'p-value']
            est.vector[i] = glmm.values['as.factor(trt2)1', 'Value']
            se.vector[i] = glmm.values['as.factor(trt2)1', 'Std.Error']
            stat.vector[i] = glmm.values['as.factor(trt2)1', 't-value']
          }
          converge.vector[i] <-
            ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE)
        }
        
        if (sigma_sq == sigma_sq2 && sigma_b_sq != sigma_b_sq2) {
          my.mod <-
            lmerTest::lmer(y ~ trt + (1 + as.factor(trt) | clust),
                           REML = FALSE,
                           data = sim.dat)
          # get the overall p-values (>Chisq)
          null.mod <-
            stats::update.formula(my.mod, y ~ 1 + (1 + as.factor(trt) |
                                                     clust))
          glmm.values = summary(my.mod)$coefficients
          pval.vector[i] = glmm.values['trt', 'Pr(>|t|)']
          est.vector[i] = glmm.values['trt', 'Estimate']
          se.vector[i] = glmm.values['trt', 'Std. Error']
          stat.vector[i] = glmm.values['trt', 't value']
          converge.vector[i] = ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == FALSE, TRUE, FALSE)
          # option to stop the function early if fits are singular
          if (poorFitOverride == FALSE) {
            if (sum(converge.vector == FALSE, na.rm = TRUE) > (nsim * .25)) {
              stop(
                "more than 25% of simulations are singular fit: check model specifications"
              )
            }
          }
        }
        
        if (sigma_sq != sigma_sq2 && sigma_b_sq == sigma_b_sq2) {
          trt2 <- unlist(trt)
          clust2 <- unlist(clust)
          oldw <- getOption("warn")
          options(warn = -1)
          my.mod <-
            try(nlme::lme(
              y ~ as.factor(trt2),
              random =  ~ 1 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ),
            silent = TRUE)
          glmm.values <-  summary(my.mod)$tTable
          # get the overall p-values (>Chisq)
          null.mod <-
            try(nlme::lme(
              y ~ 1,
              random =  ~ 1 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ),
            silent = TRUE)
          options(warn = oldw)
          pval.vector[i] = glmm.values['as.factor(trt2)1', 'p-value']
          est.vector[i] = glmm.values['as.factor(trt2)1', 'Value']
          se.vector[i] = glmm.values['as.factor(trt2)1', 'Std.Error']
          stat.vector[i] = glmm.values['as.factor(trt2)1', 't-value']
          converge.vector[i] <-
            ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE)
        }
        
        if (sigma_sq == sigma_sq2 && sigma_b_sq == sigma_b_sq2) {
          my.mod <-  lmerTest::lmer(y ~ trt + (1 | clust), REML = FALSE,
                                    data = sim.dat)
          # get the overall p-values (>Chisq)
          null.mod <- update.formula(my.mod, y ~ 1 + (1 | clust))
          glmm.values = summary(my.mod)$coefficients
          pval.vector[i] = glmm.values['trt', 'Pr(>|t|)']
          est.vector[i] = glmm.values['trt', 'Estimate']
          se.vector[i] = glmm.values['trt', 'Std. Error']
          stat.vector[i] = glmm.values['trt', 't value']
          converge.vector[i] = ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == TRUE, FALSE, TRUE)
          # option to stop the function early if fits are singular
          if (poorFitOverride == FALSE) {
            if (sum(converge.vector == FALSE, na.rm = TRUE) > (nsim * .25)) {
              stop(
                "more than 25% of simulations are singular fit: check model specifications"
              )
            }
          }
        }
      }
    }
    
    # Fit GEE (geeglm)
    # Note: there is no option for GEE with irgtt
    if (method == 'gee') {
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(y ~ trt,
                               data = sim.dat,
                               id = clust,
                               corstr = "exchangeable")
      gee.values = summary(my.mod)$coefficients
      est.vector[i] = gee.values['trt', 'Estimate']
      se.vector[i] = gee.values['trt', 'Std.err']
      stat.vector[i] = gee.values['trt', 'Wald']
      pval.vector[i] = gee.values['trt', 'Pr(>|W|)']
      converge.vector[i] <- ifelse(summary(my.mod)$error == 0, TRUE, FALSE)
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
    
    # Update simulation progress information
    if (quiet == FALSE) {
      if (i == 1) {
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
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
  }
  
  ## Output objects
  # Create object containing summary statement
  if (irgtt == FALSE) {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: Parallel Design, Continuous Outcome"
    )
  } else {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: IRGTT Design, Continuous Outcome"
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
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Arm.1' = nsubjects[1:nclusters[1]],
                       'Arm.2' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm.1" = c("n.clust" = nclusters[1]),
    "Arm.2" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame(
    'Arm.1' = c(
      'ICC' = ICC[1],
      'sigma_sq' = sigma_sq[1],
      'sigma_b_sq' = sigma_b_sq[1],
      'mu' = mu
    ),
    'Arm.2' = c(
      'ICC' = ICC2,
      'sigma_sq' = sigma_sq2,
      'sigma_b_sq' = sigma_b_sq2,
      'mu' = mu2
    )
  ))
  
  fail <- unlist(converge.vector)
  
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
      "convergence" = fail,
      "sim.data" = simulated.datasets
    ),
    class = 'crtpwr'
  )
  return(complete.output)
}