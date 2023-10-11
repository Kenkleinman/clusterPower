#' Power simulations for cluster-randomized trials: Parallel Designs, Binary Outcome
#'
#' @description
#' \loadmathjax
#'
#'
#' This function uses Monte Carlo methods (simulations) to estimate
#' power for cluster-randomized trials. Users
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#'
#' Users must specify the desired number of simulations, number of subjects per
#' cluster, number of clusters per arm, and two of the following three
#' parameters: expected probability of the outcome in one group, expected
#' probability of the outcome in the second group,
#' and expected difference in probabilities between groups.
#' Default values are provided for significance level, analytic method,
#' progress updates, and whether the simulated data sets are retained.
#'
#'
#' @param nsim Number of datasets to simulate; accepts integer. Required.
#'
#' @param nsubjects Number of subjects per cluster; accepts either a scalar
#' (implying equal cluster sizes for the two groups), a vector of length two
#' (equal cluster sizes within arm), or a vector of length \code{sum(nclusters)}
#' (unequal cluster sizes within arm). Required.
#'
#' @param nclusters Number of clusters per treatment group; accepts a single integer
#' (if there are the same number of clusters in each arm) or a vector of 2 integers
#' (if nsubjects differs between arms). If a vector of cluster sizes >2 is provided in
#' \code{nsubjects}, \code{sum(nclusters)} must match the \code{nsubjects} vector length.
#' Required.
#'
#' @param p1 Expected probability of outcome in first group.
#' @param p2 Expected probability of outcome in second group.
#'
#' @param sigma_b_sq Between-cluster variance; if sigma_b_sq2 is not specified,
#' between-cluster variances are assumed to be equal in the two arms. Accepts numeric. Required.
#' @param sigma_b_sq2 Between-cluster variance for clusters in second group. Only required if
#' between-cluster variances differ between treatment arms.
#'
#'
#' @param alpha Significance level; default = 0.05.
#' @param method Data analysis method, either generalized linear mixed effects model (GLMM)
#' or generalized estimating equations (GEE). Accepts c('glmm', 'gee'); default = 'glmm'. Required.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion
#' time, default = TRUE.
#' @param allSimData Option to include a list of all simulated datasets in the output object.
#' Default = \code{FALSE}.
#' @param nofit Option to skip model fitting and analysis and instead return a dataframe with
#' the simulated datasets. Default = \code{FALSE}.
#' @param allSimData Option to output list of all simulated datasets; default = FALSE.
#' @param nofit Option to skip model fitting and analysis and only return the simulated data.
#' Default = \code{FALSE}.
#' @param seed Option to set the seed. Default is NA.
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' @param irgtt Logical. Default = FALSE. Is the experimental design an
#' individually randomized group treatment trial? For details,
#' see ?cps.irgtt.binary.
#'
#' @return If \code{nofit = F}, a list with the following components:
#' \itemize{
#'   \item Character string indicating total number of simulations, simulation type,
#'   and number of convergent models
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (estimated statistical power),
#'   "lower.95.ci" (lower 95% confidence interval bound),
#'   "upper.95.ci" (upper 95% confidence interval bound),
#'   "Alpha" (probability of committing a Type I error or rejecting a true null),
#'   "Beta" (probability of committing a Type II error or failing to reject a false null).
#'   Note that non-convergent models are returned for review,
#'   but not included in this calculation.
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting sigma_b_sq for each group
#'   \item Vector containing user-supplied outcome probability and estimated odds ratio
#'   \item Data frame containing three estimates of ICC
#'   \item Data frame with columns:
#'   "Estimate" (Estimate of treatment effect for a given simulation),
#'   "Std.err" (Standard error for treatment effect estimate),
#'   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)),
#'   "p.value",
#'   "converge" (Did simulated model converge?)
#'   \item If allSimData = TRUE, list of data frames, each containing: "y" (Simulated response value),
#'   "trt" (Indicator for treatment group), "clust" (Indicator for cluster)
#'   \item List of warning messages produced by non-convergent models;
#'   Includes model number for cross-referencing against \code{model.estimates}
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
#'
#' @details
#'
#' The data generating model for observation \mjseqn{j} in cluster \mjseqn{i} is:
#'
#' \mjsdeqn{y_{ij} \sim \code{Bernoulli}(\frac{e^{p_1 + b_i}}{1 + e^{p_1 + b_i} }) }
#' for the first group or arm, where \mjseqn{b_i \sim N(0,\sigma_b^2)},
#' while for the second group,
#' \mjsdeqn{y_{ij} \sim \code{Bernoulli}(\frac{e^{p_2 + b_i}}{1 + e^{p_2 + b_i} }) }
#' where \mjseqn{b_i \sim N(0,\sigma_{b_2}^2)}; if
#' \mjseqn{\sigma_{b_2}^2} is not used, then the second group uses
#' \mjseqn{b_i \sim N(0,\sigma_b^2)}.
#'
#' All random terms are generated independent of one another.
#'
#'
#' Non-convergent models are not included in the calculation of exact confidence
#' intervals.
#'
#'
#' @seealso
#'
#' An intracluster correlation coefficient (ICC) for binary outcome data is
#' neither a natural parameter of the data generating model nor a function
#' of its parameters.  Several methods for calculation have been suggested
#' (Wu, Crespi, and Wong, 2012).  We provide several versions of ICCs for
#' comparison.  These can be accessed in the \code{bincalcICC()} function.
#'
#'
#'
#' @section Testing details:
#' This function has been verified against reference values from the NIH's GRT
#' Sample Size Calculator, PASS11, \code{CRTsize::n4prop}, and
#' \code{clusterPower::cpa.binary}.
#'
#' @author Alexander R. Bogdan, Alexandria C. Sakrejda
#' (\email{acbro0@@umass.edu}), and Ken Kleinman
#' (\email{ken.kleinman@@gmail.com})
#' #'
#' @references Elridge, S., Ukoumunne, O. & Carlin, J. The Intra-Cluster Correlation
#' Coefficient in Cluster Randomized Trials:
#' A Review of Definitions. International Statistical Review (2009), 77, 3, 378-394.
#' doi: 10.1111/j.1751-5823.2009.00092.x
#' @references Snjiders, T. & Bosker, R. Multilevel Analysis: an Introduction to Basic and
#' Advanced Multilevel Modelling. London, 1999: Sage.
#' @references Wu S, Crespi CM, Wong WK. Comparison of Methods for Estimating Intraclass
#' Correlation Coefficient for Binary Responses in Cancer Prevention Cluster Randomized
#' Trials. Contemp Clin Trials. 2012; 33(5): 869-880. doi:10.1016/j.cct.2012.05.004
#' London: Arnold; 2000.
#'
#' @examples
#'
#' # Estimate power for a trial with 10 clusters in each arm, 20 subjects in
#' # each cluster, with a probability of 0.8 in the first arm and 0.5 in the
#' # second arm, with a sigma_b_sq = 1 in the first arm sigma_b_sq = 1.2 in
#' # the second arm.
#'
#' \dontrun{
#' binary.sim = cps.binary(nsim = 100, nsubjects = 20,
#'   nclusters = 10, p1 = 0.8,
#'   p2 = 0.5, sigma_b_sq = 1,
#'   sigma_b_sq2 = 1.2, alpha = 0.05,
#'   method = 'glmm', allSimData = FALSE)
#' }
#'
#' # Estimate power for a trial just as above, except that in the first arm,
#' # the clusters have 10 subjects in 9 of the 10 clusters and 100 in the tenth
#' # cluster, while in the second arm all clusters have 20 subjects.
#'
#' \dontrun{
#' binary.sim2 = cps.binary(nsim = 100,
#'   nsubjects = c(c(rep(10,9),100), rep(20,10)),
#'   nclusters = 10, p1 = 0.8,
#'   p2 = 0.5, sigma_b_sq = 1,
#'   sigma_b_sq2 = 1.2, alpha = 0.05,
#'   method = 'gee', allSimData = FALSE)
#' }
#'
#'
#'
#' @export

# Define function

cps.binary = function(nsim = NULL,
                      nsubjects = NULL,
                      nclusters = NULL,
                      p1 = NULL,
                      p2 = NULL,
                      sigma_b_sq = NULL,
                      sigma_b_sq2 = NULL,
                      alpha = 0.05,
                      method = 'glmm',
                      quiet = FALSE,
                      allSimData = FALSE,
                      seed = NA,
                      nofit = FALSE,
                      poorFitOverride = FALSE,
                      lowPowerOverride = FALSE,
                      timelimitOverride = TRUE,
                      irgtt = FALSE) {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create objects to collect iteration-specific values
  est.vector <- NULL
  se.vector <- NULL
  stat.vector <- NULL
  pval.vector <- NULL
  converge.ind <- NULL
  converge.vector <- NULL
  icc2.vector <- NULL
  lmer.icc.vector <- NULL
  converge.vector <- NULL
  simulated.datasets <- list()
  warning.list <- list()
  converge <- NULL
  
  # Create progress bar
  prog.bar =  progress::progress_bar$new(
    format = "(:spin) [:bar] :percent eta :eta",
    total = nsim,
    clear = FALSE,
    width = 100
  )
  prog.bar$tick(0)
  
  # Define wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5)
    abs(x - round(x)) < tol
  
  # Define expit function
  expit = function(x)
    1 / (1 + exp(-x))
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS
  sim.data.arg.list = list(nsim, nsubjects, nclusters, sigma_b_sq)
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
  if (any(!is.wholenumber(nsubjects) | nsubjects < 1)) {
    stop(paste0("NSUBJECTS", min1.warning))
  }
  if (!is.wholenumber(nclusters) || nclusters < 1) {
    stop(paste0("NCLUSTERS", min1.warning))
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
  
  if (irgtt == FALSE) {
    # Validate sigma_b_sq, sigma_b_sq2
    min0.warning = " must be a numeric value greater than 0"
    if (!is.numeric(sigma_b_sq) || sigma_b_sq <= 0) {
      stop("sigma_b_sq", min0.warning)
    }
    if (!is.null(sigma_b_sq2) && sigma_b_sq2 <= 0) {
      stop("sigma_b_sq2", min0.warning)
    }
  }
  # Set between-cluster variances
  if (is.null(sigma_b_sq2)) {
    sigma_b_sq[2] = sigma_b_sq
  } else{
    sigma_b_sq[2] = sigma_b_sq2
  }
  
  # Validate P1, P2
  parm1.arg.list = list(p1, p2)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if (sum(parm1.args) > 1) {
    stop("Both terms must be specified: p1, p2")
  }
  
  # Validate ALPHA, METHOD, QUIET, allSimData
  if (!is.numeric(alpha) || alpha < 0) {
    stop("ALPHA", min0.warning)
  } else if (alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
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
  
  # Calculate ICC1 (sigma_b_sq / (sigma_b_sq + pi^2/3))
  icc1 = mean(sapply(1:2, function(x)
    sigma_b_sq[x] / (sigma_b_sq[x] + pi ^ 2 / 3)))
  
  # Create indicators for arm & cluster
  trt = c(rep(1, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(2, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  trt <- as.factor(trt)
  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  clust <- as.factor(clust)
  
  # Calculate log odds for each group
  logit.p1 = log(p1 / (1 - p1))
  logit.p2 = log(p2 / (1 - p2))
  
  # for irgtt option
  index <- 1
  
  ### Create simulation loop
  while (sum(converge.vector == TRUE) != nsim) {
    # Generate between-cluster effects for arm 1 and arm 2
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq[1]))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq[2]))
    
    # Create arm 1 y-value
    y0.intercept = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.0[x], length.out = nsubjects[x])))
    y0.linpred = y0.intercept + logit.p1
    y0.prob = expit(y0.linpred)
    y0 = unlist(lapply(y0.prob, function(x)
      stats::rbinom(1, 1, x)))
    if (length(table(y0)) != 2) {
      warning(print(
        "y0 is completely seperated. Repeating the random draw 1 time."
      ))
      randint.0 = stats::rnorm(nclusters[1],
                               mean = 0,
                               sd = sqrt(sigma_b_sq[1]))
      y0.intercept = unlist(lapply(1:nclusters[1], function(x)
        rep(randint.0[x], length.out = nsubjects[x])))
      y0.linpred = y0.intercept + logit.p1
      y0.prob = expit(y0.linpred)
      y0 = unlist(lapply(y0.prob, function(x)
        stats::rbinom(1, 1, x)))
    }
    
    # Create arm 2 y-value
    y1.intercept = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.linpred = y1.intercept + logit.p2
    y1.prob = expit(y1.linpred)
    y1 = unlist(lapply(y1.prob, function(x)
      stats::rbinom(1, 1, x)))
    if (length(table(y1)) != 2) {
      warning(print(
        "y1 is completely seperated. Repeating the random draw 1 time."
      ))
      randint.1 = stats::rnorm(nclusters[2],
                               mean = 0,
                               sd = sqrt(sigma_b_sq[2]))
      y1.intercept = unlist(lapply(1:nclusters[2], function(x)
        rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
      y1.linpred = y1.intercept + logit.p2
      y1.prob = expit(y1.linpred)
      y1 = unlist(lapply(y1.prob, function(x)
        stats::rbinom(1, 1, x)))
    }
    
    
    
    # Create single response vector
    y = c(y0, y1)
    
    # Create and store data frame for simulated dataset
    sim.dat = data.frame(y = y, trt = trt, clust = clust)
    if (allSimData == TRUE && nofit == FALSE) {
      simulated.datasets = append(simulated.datasets, list(sim.dat))
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
    
    # Calculate ICC2 ([P(Yij = 1, Yih = 1)] - pij * pih) / sqrt(pij(1 - pij) * pih(1 - pih))
    #icc2 = (mean(y0.prob) * mean(y1.prob) - p1*p2) / sqrt((p1 * (1 - p1)) * p2 * (1 - p2))
    icc2 = (mean(y0.prob) - p1) * (mean(y1.prob) - p2) / sqrt((p1 * (1 - p1)) * p2 * (1 - p2))
    # ^Equation above #11 (no number); Eldridge, Ukoumunne & Carlin, 2009 (p.386)
    icc2.vector = append(icc2.vector, icc2)
    
    # Calculate LMER.ICC (lmer: sigma_b_sq / (sigma_b_sq + sigma))
    if (irgtt == FALSE) {
      lmer.mod = lme4::glmer(
        y ~ trt + (1 | clust),
        data = sim.dat,
        family = stats::binomial(link = 'logit')
      )
      lmer.vcov <-
        as.numeric(as.data.frame(lme4::VarCorr(lmer.mod))[, 4:5])
      icc.val <- lmer.vcov[1] / (lmer.vcov[1] + lmer.vcov[2])
      lmer.icc.vector <- append(lmer.icc.vector, icc.val)
    }
    
    # Set warnings to OFF
    # Note: Warnings will still be stored in 'warning.list'
    options(warn = -1)
    
    start.time = Sys.time()
    
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      if (irgtt == TRUE) {
        my.mod <-
          try(MASS::glmmPQL(
            y ~ trt,
            random =  ~ 0 + trt | clust,
            data = sim.dat,
            family = stats::binomial(link = 'logit')
          ))
        if (class(my.mod) == "try-error") {
          glmm.values <- NA
          est.vector[index] <- NA
          se.vector[index] <- NA
          stat.vector[index] <- NA
          pval.vector[index] <- NA
          converge.vector[index] <- FALSE
        } else {
          glmm.values <- summary(my.mod)$tTable
          est.vector[index] <- glmm.values['trt2', 'Value']
          se.vector[index] <- glmm.values['trt2', 'Std.Error']
          stat.vector[index] <- glmm.values['trt2', 't-value']
          pval.vector[index] <- glmm.values['trt2', 'p-value']
          converge.vector[index] <- TRUE
        }
        index <- index + 1
      } else {
        my.mod = try(lme4::glmer(
          y ~ trt + (1 | clust),
          data = sim.dat,
          family = stats::binomial(link = 'logit')
        ))
        model.converge = try(my.mod)
        converge.ind = is.null(model.converge@optinfo$conv$lme4$messages)
        converge.vector = append(converge.vector, converge.ind)
        if (!isTRUE(converge.ind)) {
          model.id = paste0("Model ", length(converge.vector))
          warning.list[model.id] = list(model.converge@optinfo$conv$lme4$messages)
          glmm.values = NA
          est.vector = append(est.vector, NA)
          se.vector = append(se.vector, NA)
          stat.vector = append(stat.vector, NA)
          pval.vector = append(pval.vector, NA)
        } else {
          glmm.values = summary(my.mod)$coefficient
          est.vector = append(est.vector, glmm.values['trt2', 'Estimate'])
          se.vector = append(se.vector, glmm.values['trt2', 'Std. Error'])
          stat.vector = append(stat.vector, glmm.values['trt2', 'z value'])
          pval.vector = append(pval.vector, glmm.values['trt2', 'Pr(>|z|)'])
        }
      }
      if (poorFitOverride == FALSE &&
          length(converge.vector) > 50 &&
          sum(converge.vector == FALSE, na.rm = TRUE) > (nsim * 0.25)) {
        stop("more than 25% of simulations are singular fit: check model specifications")
      }
      
    }
    
    # Set warnings to ON
    # Note: Warnings will still be stored in 'warning.list'
    options(warn = 0)
    
    # Fit GEE (geeglm)
    if (method == "gee") {
      sim.dat = dplyr::arrange(sim.dat, clust)
      if (irgtt == FALSE) {
        my.mod = geepack::geeglm(
          y ~ trt,
          data = sim.dat,
          family = stats::binomial(link = 'logit'),
          id = clust,
          corstr = "exchangeable"
        )
      } else {
        my.mod = geepack::geeglm(
          y ~ trt + (0 + trt | clust),
          data = sim.dat,
          family = stats::binomial(link = 'logit'),
          id = clust,
          corstr = "exchangeable"
        )
      }
      gee.values = summary(my.mod)$coefficients
      est.vector = append(est.vector, gee.values['trt2', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt2', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt2', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt2', 'Pr(>|W|)'])
      converge.vector = append(converge.vector, ifelse(summary(my.mod)$error == 0, TRUE, FALSE))
    }
    
    # Print simulation start message
    
    if (length(est.vector) == 1) {
      avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      time.est = avg.iter.time * (nsim - 1) / 60
      hr.est = time.est %/% 60
      min.est = round(time.est %% 60, 3)
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
      if (min.est > 2 && timelimitOverride == FALSE) {
        stop(paste0(
          "Estimated completion time: ",
          hr.est,
          'Hr:',
          min.est,
          'Min'
        ))
      }
    }
    # Update simulation progress information
    if (quiet == FALSE) {
      prog.bar$update(sum(converge.vector == TRUE) / nsim)
      Sys.sleep(1 / 100)
    }
    # stop the loop if power is <0.5
    if (lowPowerOverride == FALSE && length(pval.vector) > 50) {
      sig.val.temp <-
        ifelse(pval.vector < alpha, 1, 0)
      pval.power.temp <-
        sum(sig.val.temp, na.rm = TRUE) / length(pval.vector)
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
  }
  # Print simulation complete message
  if (quiet == FALSE && sum(converge.vector) == nsim) {
    message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
  }
  # Governor to prevent infinite non-convergence loop
  converge.ratio <-
    sum(converge.vector == FALSE) / sum(converge.vector == TRUE)
  if (converge.ratio > 4.0 && converge.ratio != Inf) {
    stop(
      "WARNING! The number of non-convergent models exceeds the number of convergent models by a factor of 4. Consider reducing sigma_b_sq"
    )
  }
  
  ## Output objects
  # Create object containing summary statement
  if (irgtt == FALSE) {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: Simple Design, Binary Outcome. Note: ",
      sum(converge.vector == FALSE),
      " additional models were fitted to account for non-convergent simulations."
    )
  } else {
    summary.message = paste0(
      "Monte Carlo Power Estimation based on ",
      nsim,
      " Simulations: IRGTT Design, Binary Outcome. Note: Models fit using penalized quasi-likelihood."
    )
  }
  
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model',
                       gee = 'Generalized Estimating Equation')
  
  # Store model estimate output in data frame
  cps.model.est = data.frame(
    Estimate = as.vector(unlist(est.vector)),
    Std.err = as.vector(unlist(se.vector)),
    Test.statistic = as.vector(unlist(stat.vector)),
    p.value = as.vector(unlist(pval.vector)),
    converge = as.vector(unlist(converge.vector))
  )
  
  # Calculate and store power estimate & confidence intervals
  if (!is.na(any(cps.model.est$converge))) {
    cps.model.temp <- dplyr::filter(cps.model.est, converge == TRUE)
  } else {
    cps.model.temp <- cps.model.est
  }
  power.parms <- confintCalc(nsim = nsim,
                             alpha = alpha,
                             p.val = cps.model.temp[, 'p.value'])
  
  # Create object containing inputs
  p1.p2.or = round(p1 / (1 - p1) / (p2 / (1 - p2)), 3)
  p2.p1.or = round(p2 / (1 - p2) / (p1 / (1 - p1)), 3)
  inputs = t(data.frame(
    'Arm1' = c("probability" = p1, "odds.ratio" = p1.p2.or),
    'Arm2' = c("probability" = p2, 'odds.ratio' = p2.p1.or)
  ))
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Arm1' = nsubjects[1:nclusters[1]],
                       'Arm2' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm1" = c("n.clust" = nclusters[1]),
    "Arm2" = c("n.clust" = nclusters[2])
  ))
  
  if (irgtt == FALSE) {
    # Create object containing estimated ICC values
    ICC = round(t(data.frame(
      'P_h' = c('ICC' = icc1),
      'P_c' = c('ICC' = mean(icc2.vector, na.rm = TRUE)),
      'lmer' = c('ICC' = mean(lmer.icc.vector, na.rm = TRUE))
    )), 3)
    # Create object containing all ICC values
    # Note: P_h is a single calculated value. No vector to be appended.
    icc.list = data.frame('P_c' = icc2.vector,
                          'lmer' = lmer.icc.vector)
  }
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame(
    'Arm1' = c('sigma_b_sq' = sigma_b_sq[1]),
    'Arm2' = c('sigma_b_sq' = sigma_b_sq[2])
  ))
  
  # Check & governor for inclusion of simulated datasets
  # Note: If number of non-convergent models exceeds 5% of NSIM,
  # override allSimData and output all simulated data sets
  
  if (allSimData == FALSE &&
      (sum(converge.vector == FALSE) < sum(converge.vector == TRUE) * 0.05)) {
    simulated.datasets = NULL
  }
  
  # Create list containing all output (class 'crtpwr') and return
  if (irgtt == FALSE) {
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
        "inputs" = inputs,
        "ICC" = ICC,
        "icc.list" = icc.list,
        "model.estimates" = cps.model.est,
        "sim.data" = simulated.datasets,
        "warning.list" = warning.list,
        "convergence" = converge.vector
      ),
      class = "crtpwr"
    )
  } else {
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
        "inputs" = inputs,
        "model.estimates" = cps.model.est,
        "sim.data" = simulated.datasets,
        "warning.list" = warning.list,
        "convergence" = converge.vector
      ),
      class = "crtpwr"
    )
  }
  return(complete.output)
}
