#' Simulation-based power estimation for cluster-randomized trials: Parallel Designs, Count Outcome with multiple arms
#'
#'
#'
#' @description 
#' \loadmathjax
#' 
#'  
#' This function uses Monte Carlo methods (simulations) to estimate 
#' power for cluster-randomized trials for integer-valued outcomes with two or more
#' trial conditions. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, between-cluster variance, and
#' two of the following three parameters: mean event rate per unit time in one group,
#' the mean event rate per unit time in the second group, and/or the
#' mean difference in event rates between groups. Default values are provided
#' for significance level, analytic method, progress updates, and whether the simulated data sets are retained.
#'
#' Note that if all units have the same observation time, you can use the
#' mean count instead of the "mean event per unit time" in the preceding paragraph.
#'
#'
#'
#' Users must specify the desired number of simulations, number of subjects per
#' cluster, number of clusters per treatment arm, group probabilities, and the
#' between-cluster variance. Significance level, analytic method, whether progress
#' updates are displayed, poor/singular fit override, and whether or not to return the
#' simulated data may also be specified. 
#' 
#' This user-friendly function calls an internal function; the internal function 
#' can be called
#' directly by the user to return the fitted models rather than the power
#' summaries (see \code{?cps.ma.count.internal} for details).
#'
#' Users can spread the simulated data
#' generation and model fitting tasks across multiple cores using the
#' \code{cores} argument. Users should expect that parallel computing may make
#' model fitting faster than using a single core for more complicated models.
#' For simpler models, users may prefer to use single thread computing
#' (\code{cores}=1), as the processes involved in allocating memory and
#' copying data across cores also may take some time. For time-savings,
#' this function stops execution early if estimated power < 0.5 or more
#' than 25\% of models produce a singular fit or non-convergence warning
#' message, unless \code{poor.fit.override = TRUE}.
#'
#'
#'
#'
#'
#' @param nsim Number of datasets to simulate; accepts integer. Required.
#'
#' @param nsubjects Number of subjects per cluster; accepts an
#' integer (implying equal cluster sizes in all arms) if \code{narms}
#' and \code{nclusters} are provided. Alternately, a list with one integer per arm (if the
#' cluster sizes are the same within the arm), or a list of vectors where each vector represents an arm
#' and each entry in the vector is the number of subjects per cluster (if the cluster sizes are not the same
#' within the arms). Required.
#'
#' @param narms Number of trial arms; accepts integer. Required.
#'
#' @param nclusters Number of clusters per treatment group; accepts a single integer
#' (if there are the same number of clusters in each arm) or a vector of integers
#' representing the number of clusters in each arm (if nsubjects differs between arms).
#' If a list of vectors of cluster sizes is provided in \code{nsubjects}, then
#' the vector of cluster counts must match the length of the \code{nsubjects} vectors.
#' Required.
#'
#' @param counts Mean event rates per unit time for each arm; accepts a scalar
#' (if all arms have the same event rate) or
#' a vector of length \code{narms}. Required.
#'
#' @param family Distribution from which responses are simulated. Accepts Poisson
#' (\code{'poisson'}) or negative binomial (\code{'neg.binom'}); default = 'poisson'. Required.
#'
#' @param analysis Family used for data analysis; currently only applicable when \code{method = 'glmm'}.
#' Accepts c('poisson', 'neg.binom'); default = 'poisson'. Required.
#' 
#' @param negBinomSize Only used when generating simulated data from the 
#' negative binomial (family = 'neg.binom'), this is the target for number of 
#' successful trials, or the dispersion parameter (the shape parameter of the gamma 
#' mixing distribution). Must be positive and defaults to 1. Required when 
#' family = 'neg.binom'.
#' 
#' @param sigma_b_sq Between-cluster variance for each arm; accepts a scalar
#' (if all arms have the same between-cluster variance) or a vector of length
#' \code{narms}. Required.
#'
#' @param alpha The level of significance of the test, the probability of a
#' Type I error. Default = 0.05.
#'
#' @param quiet When set to \code{FALSE}, displays simulation progress and estimated
#' completion time. Default = \code{FALSE}.
#'
#' @param method Data analysis method, either generalized linear mixed effects model
#' (GLMM) or generalized estimating equations (GEE). Accepts \code{c('glmm', 'gee')};
#' default = \code{'glmm'}. Required.
#'
#' @param multi.p.method A string indicating the method to use for adjusting
#' p-values for multiple comparisons. Choose one of "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", "none". The default is
#' "bonferroni". See \code{?p.adjust} for additional details.
#'
#' @param all.sim.data Option to include a list of all simulated datasets in the output object.
#' Default = \code{FALSE}.
#'
#' @param seed Option to set the seed. Default is NULL.
#'
#' @param cores Number of cores to be used for parallel computing. Accepts a
#' string ("all"), NA (no parallel computing), or scalar value indicating
#' the number of CPUs to use. Default = NA.
#'
#' @param tdist.re Logical value indicating whether cluster-level random effects
#' should be drawn from a \mjseqn{t} distribution rather than a normal distribution.
#' Default = \code{FALSE}.
#'
#' @param poor.fit.override Option to override \code{stop()} if more than 25\% of fits fail to converge;
#' default = \code{FALSE}.
#'
#' @param low.power.override Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On \code{stop}, the power calculated from the completed simulations is printed in the
#' stop message. Default = \code{FALSE}. When \code{TRUE}, this check is ignored and the
#' calculated power is returned regardless of value.
#'
#' @param return.all.models Logical; Returns all of the fitted models, the simulated data,
#' the overall model comparisons, and the convergence report vector. This is equivalent
#' to the output of cps.ma.count.internal(). See ?cps.ma.count.internal() for details.
#'
#' @param nofit Option to skip model fitting and analysis and return the simulated data.
#' Defaults to \code{FALSE}.
#'
#' @param opt Optimizer for model fitting, from the package \code{optimx} or \code{nloptwrap}.
#' Default is 'NLOPT_LN_BOBYQA'. XX KK JN Note: This needs to be discussed more with Lexi in
#' conjuction with other optimization/singularity/convergence problems.
#'
#'
#'
#'
#' @return A list with the following components:
#' \describe{
#'   \item{power}{Data frame with columns "power" (Estimated statistical power),
#'                "lower.95.ci" (Lower 95\% confidence interval bound),
#'                "upper.95.ci" (Upper 95\% confidence interval bound).}
#'   \item{model.estimates}{Data frame with columns corresponding
#'   to each arm with descriptive suffixes as follows:
#'                   ".Estimate" (Estimate of treatment effect for a given
#'                   simulation),
#'                   "Std.Err" (Standard error for treatment effect estimate),
#'                   ".zval" (for GLMM) | ".wald" (for GEE), and
#'                   ".pval" (the p-value estimate).}
#'   \item{overall.power}{Table of F-test (when method="glmm") or chi^{2}
#'   (when method="gee") significance test results.}
#'   \item{overall.power.summary}{Summary overall power of treatment model
#'   compared to the null model.}
#'   \item{sim.data}{Produced when all.sim.data==TRUE. List of \code{nsim}
#'   data frames, each containing:
#'                   "y" (simulated response value),
#'                   "trt" (indicator for treatment group or arm), and
#'                   "clust" (indicator for cluster).}
#'   \item{model.fit.warning.percent}{Character string containing the percent
#'   of \code{nsim} in which the glmm fit was singular or failed to converge,
#'   produced only when method = "glmm" & all.sim.data = FALSE.
#'   }
#'   \item{model.fit.warning.incidence}{Vector of length \code{nsim} denoting
#'   whether or not a simulation glmm fit triggered a "singular fit" or
#'   "non-convergence" error, produced only when method = "glmm" &
#'   all.sim.data=TRUE.
#'   }
#'   }
#'   If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#'
#'
#' @details
#'
#' If \code{family = 'poisson'}, the data generating model is:
#' \mjsdeqn{y_{ijk} \sim \text{Poisson}(e^{c_k + b_{jk}}) }
#' for observation \mjseqn{i}, in cluster \mjseqn{j}, in treatment arm \mjseqn{k}, where \mjseqn{b_{jk}\sim N(0,\sigma^2_{b_{k}})}.
#'
#' If \code{family = 'neg.bin'}, the data generating model, using the
#' alternative parameterization of the negative binomial distribution
#' detailed in \code{stats::rnbinom}, is:
#'
#' \mjsdeqn{y_{ijk} \sim \text{NB}(\mu = e^{c_k + b_{jk}}, \text{size} = 1) }
#'
#' for observation \mjseqn{i}, in cluster \mjseqn{j}, in treatment arm \mjseqn{k}, where \mjseqn{b_{jk}\sim N(0,\sigma^2_{b_{k}})}.
#'
#'
#' Non-convergent models are not included in the calculation of exact confidence
#' intervals.
#'
#' For complicated models, we recommend using parallel processing with the \code{cores="all"} argument.
#' For simpler models, users may prefer to use single thread computing
#' (\code{cores}=1), as the processes involved in allocating memory and
#' copying data across cores also may take some time.
#'
#' By default, this function stops execution early if estimated power < 0.5 or if more
#' than 25\% of models produce a singular fit or non-convergence warning. In some cases, users
#' may want to ignore singularity warnings (see \code{?isSingular}) by setting \code{poor.fit.override = TRUE}.
#'
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Alexander R. Bogdan
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @examples
#' 
#' # For a 3-arm trial with 4, 4, and 5 clusters in each arm, respectively, 
#' # specify the number of subjects in each cluster with 3 vectors in a list, 
#' # each vector representing a study arm. For each cluster, in no particular
#' # order, denote the number of subjects. In this example, the first arm 
#' # contains 150, 200, 50, and 100 subjects in each of the 4 clusters. The second
#' # arm contains 50, 150, 210, and 100 subjects in each of 4 clusters, while 
#' # the third arm contains 70, 200, 150, 50, and 100 subjects in each of 5 
#' # clusters. The expected outcomes for each arm are 10, 55, and 65, and 
#' # the sigma_b_sq values are 1, 1, and 2, respectively. Assuming 
#' # seed = 123, the overall power for this trial should be 0.81.
#' 
#' \dontrun{
#' nsubjects.example <- list(c(150, 200, 50, 100), c(50, 150, 210, 100), c(70, 200, 150, 50, 100))
#' counts.example <- c(10, 55, 65)
#' sigma_b_sq.example <- c(1, 1, 2)
#'
#' count.ma.rct.unbal <- cps.ma.count(nsim = 100,
#'                             nsubjects = nsubjects.example,
#'                             counts = counts.example,
#'                             sigma_b_sq = sigma_b_sq.example,
#'                             alpha = 0.05, seed = 123)
#'}
#'
#' # For a different trial with 4 arms, each arm has 4 clusters which 
#' # each contain 100 subjects. Expected counts for each arm are 30 
#' # for the first arm, 35 for the second, 70 for the third, and 40
#' # for the fourth. Similarly, sigma_b_sq for each arm are 1 
#' # for the first arm, 1.2 for the second, 1 for the third, and 0.9
#' # for the fourth. Assuming seed = 123, the overall power for this 
#' # trial should be 0.84
#'
#' \dontrun{
#' count.ma.rct.bal <- cps.ma.count(nsim = 10, nsubjects = 100, narms = 4,
#'                             nclusters = 25, counts = c(30, 35, 70, 40),
#'                             sigma_b_sq = c(1, 1.2, 1, 0.9), seed = 123)
#'}
#' @export

cps.ma.count <- function(nsim = 1000,
                         nsubjects = NULL,
                         narms = NULL,
                         nclusters = NULL,
                         counts = NULL,
                         family = "poisson",
                         analysis = "poisson",
                         negBinomSize = 1,
                         sigma_b_sq = NULL,
                         alpha = 0.05,
                         quiet = FALSE,
                         method = 'glmm',
                         multi.p.method = "bonferroni",
                         all.sim.data = FALSE,
                         seed = NA,
                         cores = NA,
                         tdist.re = FALSE,
                         poor.fit.override = FALSE,
                         low.power.override = FALSE,
                         return.all.models = FALSE,
                         nofit = FALSE,
                         opt = "NLOPT_LN_BOBYQA") {
  # use this later to determine total elapsed time
  start.time <- Sys.time()
  
  # create narms and nclusters if not provided directly by user
  if (isTRUE(is.list(nsubjects))) {
    # create narms and nclusters if not supplied by the user
    if (is.null(narms)) {
      narms <- length(nsubjects)
    }
    if (is.null(nclusters)) {
      nclusters <- vapply(nsubjects, length, 0)
    }
  }
  if (length(nclusters) == 1 & !isTRUE(is.list(nsubjects))) {
    nclusters <- rep(nclusters, narms)
  }
  if (length(nclusters) > 1 & length(nsubjects) == 1) {
    narms <- length(nclusters)
  }
  
  # input validation steps
  if (!is.wholenumber(nsim) || nsim < 1 || length(nsim) > 1) {
    stop("nsim must be a positive integer of length 1.")
  }
  if (isTRUE(is.null(nsubjects))) {
    stop("nsubjects must be specified. See ?cps.ma.count for help.")
  }
  if (length(nsubjects) == 1 & !isTRUE(is.numeric(nclusters))) {
    stop("When nsubjects is scalar, user must supply nclusters (clusters per arm)")
  }
  if (length(nsubjects) == 1 & length(nclusters) == 1 &
      !isTRUE(is.numeric(narms))) {
    stop("User must provide narms when nsubjects and nclusters are both scalar.")
  }
  
  # nclusters must be whole numbers
  if (sum(is.wholenumber(nclusters) == FALSE) != 0 ||
      nclusters < 1) {
    stop("nclusters must be postive integer values.")
  }
  
  # nsubjects must be whole numbers
  if (sum(is.wholenumber(unlist(nsubjects)) == FALSE) != 0 ||
      unlist(nsubjects) < 1) {
    stop("nsubjects must be positive integer values.")
  }
  
  # Create nsubjects structure from narms and nclusters when nsubjects is scalar
  if (length(nsubjects) == 1) {
    str.nsubjects <- lapply(nclusters, function(x)
      rep(nsubjects, x))
  } else {
    str.nsubjects <- nsubjects
  }
  
  # allows for counts, sigma_b_sq to be entered as scalar
  if (length(sigma_b_sq) == 1) {
    sigma_b_sq <- rep(sigma_b_sq, narms)
  }
  if (length(counts) == 1) {
    counts <- rep(counts, narms)
  }
  
  if (length(counts) != narms) {
    stop(
      "Length of counts must equal narms, or be provided as a scalar if counts for all arms are equal."
    )
  }
  
  if (length(sigma_b_sq) != narms) {
    stop(
      "Length of variance parameters sigma_b_sq must equal narms, or be provided as a scalar
         if sigma_b_sq for all arms are equal."
    )
  }
  
  if (narms < 3) {
    message("Warning: LRT significance not calculable when narms<3. Use cps.count() instead.")
  }
  
  # validateVariance(dist="bin", alpha=alpha, ICC=NA, sigma=NA,
  #                   sigma_b=sigma_b_sq, ICC2=NA, sigma2=NA,
  #                   sigma_b2=NA, method=method, quiet=quiet,
  #                   all.sim.data=all.sim.data,
  #                   poor.fit.override=poor.fit.override,
  #                   cores=cores)
  
  # Set warnings to OFF
  # Note: Warnings will still be stored in 'warning.list'
  options(warn = -1)
  
  # run the simulations
  count.ma.rct <- cps.ma.count.internal(
    nsim = nsim,
    str.nsubjects = str.nsubjects,
    counts = counts,
    sigma_b_sq = sigma_b_sq,
    alpha = alpha,
    quiet = quiet,
    method = method,
    all.sim.data = all.sim.data,
    seed = seed,
    poor.fit.override = poor.fit.override,
    low.power.override = low.power.override,
    tdist = tdist.re,
    cores = cores,
    family = family,
    analysis = analysis,
    negBinomSize = negBinomSize,
    nofit = nofit,
    opt = opt
  )
  
  # Set warnings to ON
  # Note: Warnings will still be stored in 'warning.list'
  options(warn = 0)
  
  #option to return simulated data only
  if (nofit == TRUE || return.all.models == TRUE) {
    return(count.ma.rct)
  }
  
  # Create object containing summary statement
  summary.message = paste0(
    "Monte Carlo Power Estimation based on ",
    nsim,
    " Simulations: Parallel Design, Poisson Outcome, ",
    narms,
    " Arms."
  )
  
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model',
                       gee = 'Generalized Estimating Equation')
  
  # Create object containing group-specific variance parameters
  armnames <- vector(mode = "character", length = narms)
  for (i in 1:narms) {
    armnames[i] <- paste0("Arm.", i)
  }
  
  var.parms = data.frame('sigma_b_sq' = sigma_b_sq,
                         "counts" = counts)
  rownames(var.parms) <- armnames
  # Create object containing group-specific cluster sizes
  names(str.nsubjects) <- armnames
  cluster.sizes <- str.nsubjects
  
  models <- count.ma.rct[[1]]
  
  #Organize output for GLMM
  if (method == "glmm") {
    Estimates = matrix(NA, nrow = nsim, ncol = narms)
    std.error = matrix(NA, nrow = nsim, ncol = narms)
    z.val = matrix(NA, nrow = nsim, ncol = narms)
    p.val = matrix(NA, nrow = nsim, ncol = narms)
    
    for (i in 1:nsim) {
      Estimates[i,] <- models[[i]][[10]][, 1]
      std.error[i,] <- models[[i]][[10]][, 2]
      z.val[i,] <- models[[i]][[10]][, 3]
      p.val[i,] <-
        p.adjust(models[[i]][[10]][, 4], method = multi.p.method)
    }
    
    # Organize the row/col names for the model estimates output
    keep.names <- rownames(models[[1]][[10]])
    
    names.Est <- rep(NA, narms)
    names.st.err <- rep(NA, narms)
    names.zval <- rep(NA, narms)
    names.pval <- rep(NA, narms)
    names.power <- rep(NA, narms)
    
    for (i in 1:length(keep.names)) {
      names.Est[i] <- paste(keep.names[i], ".Estimate", sep = "")
      names.st.err[i] <- paste(keep.names[i], ".Std.Err", sep = "")
      names.zval[i] <- paste(keep.names[i], ".zval", sep = "")
      names.pval[i] <- paste(keep.names[i], ".pval", sep = "")
      names.power[i] <- paste(keep.names[i], ".power", sep = "")
    }
    colnames(Estimates) <- names.Est
    colnames(std.error) <- names.st.err
    colnames(z.val) <- names.zval
    colnames(p.val) <- names.pval
    
    if (narms > 2) {
      # Organize the LRT output
      LRT.holder <-
        matrix(
          unlist(count.ma.rct[[2]]),
          ncol = 3,
          nrow = nsim,
          byrow = TRUE,
          dimnames = list(seq(1:nsim),
                          colnames(count.ma.rct[[2]][[1]]))
        )
      
      # Proportion of times P(>F)
      sig.LRT <-  ifelse(LRT.holder[, 3] < alpha, 1, 0)
      LRT.holder.abbrev <- sum(sig.LRT)
    }
    
    # Calculate and store power estimate & confidence intervals
    sig.val <-  ifelse(p.val < alpha, 1, 0)
    pval.power <- apply(sig.val, 2, sum)
    
    converged <- as.vector(rep(NA, times = nsim))
    for (i in 1:nsim) {
      converged[i] <-
        ifelse(is.null(count.ma.rct[[1]][[i]]$optinfo$conv$lme4$messages),
               TRUE,
               FALSE)
    }
    cps.model.temp <- data.frame(converged, p.val)
    colnames(cps.model.temp)[1] <- "converged"
    cps.model.temp2 <-
      dplyr::filter(cps.model.temp, converged == TRUE)
    if (isTRUE(nrow(cps.model.temp2) < (.25 * nsim))) {
      warning(paste0(
        nrow(cps.model.temp2),
        " models converged. Check model parameters."
      ),
      immediate. = TRUE)
    }

    # Calculate and store power estimate & confidence intervals
    power.parms <- confintCalc(alpha = alpha,
                                p.val = as.vector(cps.model.temp2[, 3:length(cps.model.temp2)]))
    
    # Store simulation output in data frame
    ma.model.est <-
      data.frame(Estimates, std.error, z.val, p.val, count.ma.rct[[3]])
    ma.model.est <-
      ma.model.est[, -grep('.*ntercept.*', names(ma.model.est))]
    
    # performance messages
    total.est <-
      as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
    hr.est <-  total.est %/% 3600
    min.est <-  total.est %/% 60
    sec.est <-  round(total.est %% 60, 0)
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
    
    ## Output objects for GLMM
    # Create list containing all output (class 'crtpwr.ma') and return
    if (all.sim.data == TRUE && return.all.models == FALSE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est,
          "convergence" = count.ma.rct[[3]],
          "sim.data" = count.ma.rct[[4]]
        ),
        class = 'crtpwr.ma'
      )
    }
    
    if (return.all.models == TRUE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est,
          "convergence" = count.ma.rct[[3]],
          "sim.data" = count.ma.rct[[4]],
          "all.models" <-  count.ma.rct
        ),
        class = 'crtpwr.ma'
      )
    }
    if (return.all.models == FALSE && all.sim.data == FALSE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est,
          "convergence" = count.ma.rct[[3]]
        ),
        class = 'crtpwr.ma'
      )
    }
  } # end of GLMM options
  
  #Organize output for GEE method
  if (method == "gee") {
    # Organize the output
    Estimates = matrix(NA, nrow = nsim, ncol = narms)
    std.error = matrix(NA, nrow = nsim, ncol = narms)
    Wald = matrix(NA, nrow = nsim, ncol = narms)
    Pr = matrix(NA, nrow = nsim, ncol = narms)
    
    for (i in 1:nsim) {
      Estimates[i,] <- models[[i]]$coefficients[, 1]
      std.error[i,] <- models[[i]]$coefficients[, 2]
      Wald[i,] <- models[[i]]$coefficients[, 3]
      Pr[i,] <- models[[i]]$coefficients[, 4]
    }
    
    # Organize the row/col names for the output
    keep.names <- rownames(models[[1]]$coefficients)
    
    names.Est <- rep(NA, length(narms))
    names.st.err <- rep(NA, length(narms))
    names.wald <- rep(NA, length(narms))
    names.pval <- rep(NA, length(narms))
    names.power <- rep(NA, length(narms))
    
    for (i in 1:length(keep.names)) {
      names.Est[i] <- paste(keep.names[i], ".Estimate", sep = "")
      names.st.err[i] <- paste(keep.names[i], ".Std.Err", sep = "")
      names.wald[i] <- paste(keep.names[i], ".wald", sep = "")
      names.pval[i] <- paste(keep.names[i], ".pval", sep = "")
      names.power[i] <- paste(keep.names[i], ".power", sep = "")
    }
    colnames(Estimates) <- names.Est
    colnames(std.error) <- names.st.err
    colnames(Wald) <- names.wald
    colnames(Pr) <- names.pval
    
    # Organize the LRT output
    LRT.holder <-
      matrix(
        unlist(count.ma.rct[[2]]),
        ncol = 3,
        nrow = nsim,
        byrow = TRUE,
        dimnames = list(seq(1:nsim),
                        c("Df", "X2", "P(>|Chi|)"))
      )
    
    # Proportion of times P(>F)
    sig.LRT <-  ifelse(LRT.holder[, 3] < alpha, 1, 0)
    LRT.holder.abbrev <- sum(sig.LRT) / nsim
    
    # Calculate and store power estimate & confidence intervals
    power.parms <- confint.calc(nsim = nsim,
                                alpha = alpha,
                                p.val = Pr[, 2:narms])
    
    
    # Store GEE simulation output in data frame
    ma.model.est <-  data.frame(Estimates, std.error, Wald, Pr)
    ma.model.est <-
      ma.model.est[, -grep('.*ntercept.*', names(ma.model.est))]
    
    # performance messages
    total.est <-
      as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
    hr.est <-  total.est %/% 3600
    min.est <-  total.est %/% 60
    sec.est <-  round(total.est %% 60, 0)
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
    
    # Create list containing all output (class 'crtpwr.ma') and return
    if (all.sim.data == TRUE & return.all.models == FALSE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est,
          "sim.data" = count.ma.rct[[3]]
        ),
        class = 'crtpwr.ma'
      )
    }
    if (return.all.models == TRUE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est,
          "all.models" <-  count.ma.rct
        ),
        class = 'crtpwr.ma'
      )
    }
    if (return.all.models == FALSE && all.sim.data == FALSE) {
      complete.output = structure(
        list(
          "overview" = summary.message,
          "nsim" = nsim,
          "power" =  power.parms,
          "beta" = power.parms['Beta'],
          "overall.power2" =
            prop_H0_rejection(
              alpha = alpha,
              nsim = nsim,
              LRT.holder.abbrev = LRT.holder.abbrev
            ),
          "overall.power" = LRT.holder,
          "method" = long.method,
          "alpha" = alpha,
          "cluster.sizes" = cluster.sizes,
          "n.clusters" = nclusters,
          "variance.parms" = var.parms,
          "counts" = counts,
          "model.estimates" = ma.model.est
        ),
        class = 'crtpwr.ma'
      )
    }
  }
  return(complete.output)
}# end of fxn