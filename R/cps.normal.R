#' Power simulations for cluster-randomized trials: Simple Designs, Continuous Outcome.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatments, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' Non-convergent models are not included in the calculation of exact confidence 
#' intervals.
#' 
#' @section Testing details:   
#' This function has been verified against reference values from the NIH's GRT 
#' Sample Size Calculator, PASS11, \code{CRTsize::n4means}, and 
#' \code{clusterPower::cpa.normal}.
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (equal cluster sizes, both groups), 
#' a vector of length two (equal cluster sizes within groups), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within groups) (required).
#' @param nclusters Number of clusters per group; accepts single integer or vector of length 2 for unequal number 
#' of clusters per treatment group (required)
#' @param difference Expected absolute treatment effect; accepts numeric (required).
#' At least 2 of the following must be specified:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma_sq Within-cluster variance; accepts numeric
#' @param sigma_b_sq Between-cluster variance; accepts numeric
#' 
#' If clusters differ between treatment groups, at least 2 of the following 
#' must be specified:
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma_sq2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b_sq2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' @param seed Option to set the seed. Default is NA.
#' @param irgtt Logical. Is the experimental design an individually randomized 
#' group treatment trial? For details, see ?cps.irgtt.normal.
#' @param poor.fit.override Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item Character string indicating total number of simulations and simulation type
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound).
#'                Note that non-convergent models are returned for review, 
#'                but not included in this calculation.
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters in each treatment group
#'   \item Data frame reporting ICC for Treatment/Non-Treatment groups
#'   \item Vector containing expected difference between groups based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)
#'                   }
#' 
#' 
#' @examples 
#' \dontrun{
#' normal.sim = cps.normal(nsim = 1000, nsubjects = 50, nclusters = 10, difference = 3.75,
#'                         ICC = 0.3, sigma_sq = 20,
#'                         alpha = 0.05, method = 'glmm', 
#'                         quiet = FALSE, all.sim.data = FALSE)
#' }
#' 
#' @author Alexander R. Bogdan, Alexandria C. Sakrejda 
#' (\email{acbro0@@umass.edu}), and Ken Kleinman 
#' (\email{ken.kleinman@@gmail.com})
#' 
#' @export



cps.normal = function(nsim = NULL,
                      nsubjects = NULL,
                      nclusters = NULL,
                      difference = NULL,
                      ICC = NULL,
                      sigma_sq = NULL,
                      sigma_b_sq = NULL,
                      ICC2 = NULL,
                      sigma_sq2 = NULL,
                      sigma_b_sq2 = NULL,
                      alpha = 0.05,
                      method = 'glmm',
                      quiet = FALSE,
                      all.sim.data = FALSE,
                      seed = NA,
                      poor.fit.override = FALSE,
                      irgtt = FALSE) {
  # option for reproducibility
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create vectors to collect iteration-specific values
  est.vector = NULL
  se.vector = NULL
  stat.vector = NULL
  pval.vector = NULL
  # This container keeps track of how many models failed to converge
  converge.vector <- NULL
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
  
  # Validate NSIM, NCLUSTERS, NSUBJECTS
  sim.data.arg.list = list(nsim, nclusters, nsubjects, difference)
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
  if (length(nclusters) > 2) {
    stop(
      "NCLUSTERS can only be a scalar (equal # of clusters per group) or a vector of length 2 (unequal # of clusters per group)"
    )
  }
  
  # Set cluster sizes for treatment arm (if not already specified)
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
      length(nsubjects) != 1 && length(nsubjects) != sum(nclusters)) {
    stop(
      "A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS"
    )
  }

  ## Create variance parameters
  # sigma_b_sq, sigma_sq, ICC
  if (!is.null(c(ICC, sigma_sq)) && is.null(sigma_b_sq)) {
    sigma_b_sq = ICC * sigma_sq / (1 - ICC)
  }
  if (!is.null(c(ICC, sigma_b_sq)) && is.null(sigma_sq)) {
    sigma_sq = sigma_b_sq / ICC - sigma_b_sq
  }
  if (!is.null(c(sigma_sq, sigma_b_sq)) && is.null(ICC)) {
    ICC = sigma_b_sq / (sigma_b_sq + sigma_sq)
  }
  # sigma_b_sq2, sigma_sq2, ICC2
  if (!is.null(c(ICC2, sigma_sq2)) && is.null(sigma_b_sq2)) {
    sigma_b_sq2 = ICC2 * sigma_sq2 / (1 - ICC2)
  }
  if (!is.null(c(ICC2, sigma_b_sq2)) && is.null(sigma_sq2)) {
    sigma_sq2 = sigma_b_sq2 / ICC2 - sigma_b_sq2
  }
  if (!is.null(c(sigma_sq2, sigma_b_sq2)) && is.null(ICC2)) {
    ICC2 = sigma_b_sq2 / (sigma_b_sq2 + sigma_sq2)
  }
  
  # Set within/between cluster variances & ICC for treatment group (if not already specified)
  if (isTRUE(is.null(sigma_sq2))) {
    sigma_sq2 <- sigma_sq
  }
  if (isTRUE(is.null(sigma_b_sq2))) {
    sigma_b_sq2 <- sigma_b_sq
  }
  if (isTRUE(is.null(ICC2))) {
    ICC2 <- ICC
  }
  
  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value greater than 0"
  if (!is.numeric(difference) || difference < 0) {
    stop("DIFFERENCE", min0.warning)
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate ICC, sigma_sq, sigma_b_sq, ICC2, sigma_sq2, sigma_b_sq2
  
  parm1.arg.list = list(ICC, sigma_sq, sigma_b_sq)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if (sum(parm1.args) > 1) {
    stop("At least two of the following terms must be specified: ICC, sigma_sq, sigma_b_sq")
  }
  if (round(ICC, 2) != round((sigma_b_sq / (sigma_b_sq + sigma_sq)), 2)) {
    stop("At least one of the following terms has been misspecified: ICC, sigma_sq, sigma_b_sq")
  }
  
  parm2.arg.list = list(ICC2, sigma_sq2, sigma_b_sq2)
  parm2.args = unlist(lapply(parm2.arg.list, is.null))
  if (sum(parm2.args) > 1 && sum(parm2.args) != 3) {
    stop(
      "At least two of the following terms must be provided to simulate treatment-specific
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
  if (!is.logical(all.sim.data)) {
    stop(
      "ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output"
    )
  }
  
  # Create indicators for treatment group & cluster
  trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(1, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  
  # Create simulation loop
  for (i in 1:nsim) {
    # Generate between-cluster effects for non-treatment and treatment
    randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq))
    randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq2))
    
    # Create non-treatment y-value
    y0.bclust = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.0[x], length.out = nsubjects[x])))
    y0.wclust = unlist(lapply(nsubjects[1:nclusters[1]], function(x)
      stats::rnorm(
        x, mean = 0, sd = sqrt(sigma_sq)
      )))
    y.0 = y0.bclust + y0.wclust
    
    # Create treatment y-value
    y1.bclust = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.wclust = unlist(lapply(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])],
                              function(x)
                                stats::rnorm(x, mean = difference, sd = sqrt(sigma_sq2))))
    y.1 = y1.bclust + y1.wclust
    
    # Create single response vector
    y = c(y.0, y.1)
    
    # Create data frame for simulated dataset
    sim.dat = data.frame(y = y, trt = trt, clust = clust)
    if (all.sim.data == TRUE) {
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # trt and clust are re-coded as trt2 and clust2 to work nicely with lme. This can be changed later.
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
            ))
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
            ))
          pval.vector = append(pval.vector, glmm.values['as.factor(trt2)1', 'p-value'])
          est.vector = append(est.vector, glmm.values['as.factor(trt2)1', 'Value'])
          se.vector = append(se.vector, glmm.values['as.factor(trt2)1', 'Std.Error'])
          stat.vector = append(stat.vector, glmm.values['as.factor(trt2)1', 't-value'])
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
            stats::update.formula(my.mod, y ~ 1 + (0 + as.factor(trt) | clust))
          glmm.values = summary(my.mod)$coefficients
          pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|t|)'])
          est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
          se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
          stat.vector = append(stat.vector, glmm.values['trt', 't value'])
          converge.vector = append(converge.vector, ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == FALSE, TRUE))
          # option to stop the function early if fits are singular
          if (poor.fit.override == FALSE) {
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
          my.mod <-
            try(nlme::lme(
              y ~ as.factor(trt2),
              random =  ~ 1 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ))
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
            ))
          pval.vector = append(pval.vector, glmm.values['as.factor(trt2)1', 'p-value'])
          est.vector = append(est.vector, glmm.values['as.factor(trt2)1', 'Value'])
          se.vector = append(se.vector, glmm.values['as.factor(trt2)1', 'Std.Error'])
          stat.vector = append(stat.vector, glmm.values['as.factor(trt2)1', 't-value'])
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
            stats::update.formula(my.mod, y ~ 1 + (1 + as.factor(trt) | clust))
          glmm.values = summary(my.mod)$coefficients
          pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|t|)'])
          est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
          se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
          stat.vector = append(stat.vector, glmm.values['trt', 't value'])
          converge.vector = append(converge.vector, ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == FALSE, TRUE))
          # option to stop the function early if fits are singular
          if (poor.fit.override == FALSE) {
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
          my.mod <-
            try(nlme::lme(
              y ~ as.factor(trt2),
              random =  ~ 1 + as.factor(trt2) | clust2,
              weights = nlme::varIdent(form =  ~ 1 |
                                         as.factor(trt2)),
              method = "ML",
              control = nlme::lmeControl(opt = 'optim')
            ))
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
            ))
          pval.vector = append(pval.vector, glmm.values['as.factor(trt2)1', 'p-value'])
          est.vector = append(est.vector, glmm.values['as.factor(trt2)1', 'Value'])
          se.vector = append(se.vector, glmm.values['as.factor(trt2)1', 'Std.Error'])
          stat.vector = append(stat.vector, glmm.values['as.factor(trt2)1', 't-value'])
          converge.vector[i] <-
            ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE)
        }
        
        if (sigma_sq == sigma_sq2 && sigma_b_sq == sigma_b_sq2) {
          my.mod <-  lmerTest::lmer(y ~ trt + (1 | clust), REML = FALSE,
                                    data = sim.dat)
          # get the overall p-values (>Chisq)
          null.mod <- update.formula(my.mod, y ~ 1 + (1 | clust))
          glmm.values = summary(my.mod)$coefficients
          pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|t|)'])
          est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
          se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
          stat.vector = append(stat.vector, glmm.values['trt', 't value'])
          converge.vector = append(converge.vector, ifelse(any(
            grepl("singular",
                  my.mod@optinfo$conv$lme4$messages)
          ) == TRUE, FALSE, TRUE))
          # option to stop the function early if fits are singular
          if (poor.fit.override == FALSE) {
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
      est.vector = append(est.vector, gee.values['trt', 'Estimate'])
      se.vector = append(se.vector, gee.values['trt', 'Std.err'])
      stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
      pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
    }
    
    # Update simulation progress information
    if (quiet == FALSE) {
      if (i == 1) {
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 0)
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
      " Simulations: Simple Design, Continuous Outcome"
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
  power.parms <- confint.calc(
    nsim = nsim,
    alpha = alpha,
    p.val = cps.model.temp[, 'p.value'],
    names.power = c("trt")
  )
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Non.Treatment' = nsubjects[1:nclusters[1]],
                       'Treatment' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Non.Treatment" = c("n.clust" = nclusters[1]),
    "Treatment" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing group-specific variance parameters
  var.parms = t(data.frame(
    'Non.Treatment' = c(
      'ICC' = ICC[1],
      'sigma_sq' = sigma_sq[1],
      'sigma_b_sq' = sigma_b_sq[1]
    ),
    'Treatment' = c(
      'ICC' = ICC2,
      'sigma_sq' = sigma_sq2,
      'sigma_b_sq' = sigma_b_sq2
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
      "inputs" = difference,
      "model.estimates" = cps.model.est,
      "convergence" = fail,
      "sim.data" = simulated.datasets
    ),
    class = 'crtpwr'
  )
  return(complete.output)
}
