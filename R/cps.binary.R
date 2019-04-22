#' Power simulations for cluster-randomized trials: Simple Designs, Binary Outcome.
#'
#' This function utilizes iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation for binary outcomes.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, two of the following three terms: 
#' expected probability of outcome in non-treatment group, expected probability of 
#' outcome in treatment group, expected difference in probabilities between groups
#' ; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' The following equations are used to estimate intra-cluster correltation coefficients:
#' P_h: \deqn{ICC = \frac{\sigma_{b}}{\sigma_{b} + \pi^{2}/3}}
#' P_c: \deqn{ICC = \frac{P(Y_{ij} = 1, Y_{ih} = 1) - \pi_{j}\pi_{h}}{\sqrt{\pi_{j}(1 - \pi_{j})\pi_{h}(1 - \pi_{h})}}}
#' P_lmer: \deqn{ICC = \frac{\sigma_{b}}{\sigma_{b} + \sigma_{w}}}
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters per treatment group; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified when using expected probabilities:
#' @param p1 Expected probability of outcome in non-treatment group
#' @param p2 Expected probability of outcome in treatment group
#' @param p.diff Expected difference in probability of outcome between groups, defined as p.diff = p1 - p2
#' At least 2 of the following 3 arguments must be specified when using expected odds ratios:
#' @param or1 Expected odds ratio for outcome in non-treatment group
#' @param or2 Expected odds ratio for outcome in treatment group
#' @param or.diff Expected difference in odds ratio for outcome between groups, defined as or.diff = or1 - or2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal for both groups. Accepts numeric.
#' If between cluster variances differ between treatment groups, sigma_b2 must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time, default is TRUE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE
#' @param seed Option to set the seed. Default is NA.
#'  
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, simulation type, and number of convergent models
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'   "lower.95.ci" (Lower 95% confidence interval bound), "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting sigma_b for each group
#'   \item Vector containing expected difference in probabilities based on user inputs
#'   \item Data frame containing three estimates of ICC
#'   \item Data frame with columns: "Estimate" (Estimate of treatment effect for a given simulation), 
#'   "Std.err" (Standard error for treatment effect estimate), "Test.statistic" (z-value (for GLMM) or 
#'   Wald statistic (for GEE)), "p.value", "converge" (Did simulated model converge?), 
#'   "sig.val" (Is p-value less than alpha?)
#'   \item List of data frames, each containing: "y" (Simulated response value), 
#'   "trt" (Indicator for treatment group), "clust" (Indicator for cluster)
#'   \item List of warning messages produced by non-convergent models; 
#'   Includes model number for cross-referencing against \code{model.estimates}
#' }
#' 
#' @author Alexander R. Bogdan, Alexandria C. Sakrejda, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @references Snjiders, T. & Bosker, R. Multilevel Analysis: an Introduction to Basic and Advanced Multilevel Modelling. London, 1999: Sage.
#' @references Elridge, S., Ukoumunne, O. & Carlin, J. The Intra-Cluster Correlation Coefficient in Cluster Randomized Trials: 
#' A Review of Definitions. International Statistical Review (2009), 77, 3, 378-394. doi: 10.1111/j.1751-5823.2009.00092.x
#' 
#' @examples 
#' \dontrun{
#' binary.sim = cps.binary(nsim = 100, nsubjects = 20, nclusters = 10, p1 = 0.5,
#'                         p2 = 0.2, sigma_b = 1, sigma_b2 = 1, alpha = 0.05, 
#'                         method = 'glmm', all.sim.data = FALSE)
#' }
#'
#' @export

# Define function
cps.binary = function(nsim = NULL, nsubjects = NULL, nclusters = NULL, p.diff = NULL,
                        p1 = NULL, p2 = NULL, or1 = NULL, or2 = NULL, or.diff = NULL, 
                        sigma_b = NULL, sigma_b2 = NULL, alpha = 0.05, method = 'glmm', 
                      quiet = TRUE, all.sim.data = FALSE, seed = NA, irgtt = FALSE){
  if (!is.na(seed)){
  set.seed(seed = seed)
  }
    # Create objects to collect iteration-specific values
    est.vector = NULL
    se.vector = NULL
    stat.vector = NULL
    pval.vector = NULL
    converge.ind = NULL
    converge.vector = NULL
    icc2.vector = NULL
    lmer.icc.vector = NULL
    simulated.datasets = list()
    warning.list = list()
    start.time = Sys.time()
    
    # Create progress bar
    prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta", 
                                 total = nsim, clear = FALSE, width = 100)
    prog.bar$tick(0)
    
    # Define wholenumber function
    is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    # Define expit function
    expit = function(x)  1 / (1 + exp(-x))
    
    # Validate NSIM, NSUBJECTS, NCLUSTERS
    sim.data.arg.list = list(nsim, nsubjects, nclusters, sigma_b)
    sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
    if(sum(sim.data.args) > 0){
      stop("NSIM, NSUBJECTS, NCLUSTERS & SIGMA_B must all be specified. Please review your input values.")
    }
    min1.warning = " must be an integer greater than or equal to 1"
    if(!is.wholenumber(nsim) || nsim < 1){
      stop(paste0("NSIM", min1.warning))
    }
    if(!is.wholenumber(nsubjects) || nsubjects < 1){
      stop(paste0("NSUBJECTS", min1.warning))
    }
    if(!is.wholenumber(nclusters) || nclusters < 1){
      stop(paste0("NCLUSTERS", min1.warning))
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
    
    # Validate SIGMA_B, SIGMA_B2
    min0.warning = " must be a numeric value greater than 0"
    if(!is.numeric(sigma_b) || sigma_b <= 0){
      stop("SIGMA_B", min0.warning)
    }
    if(!is.null(sigma_b2) && sigma_b2 <= 0){
      stop("SIGMA_B2", min0.warning)
    }
    # Set between-cluster variances
    if(is.null(sigma_b2)){
      sigma_b[2] = sigma_b
    }else{
      sigma_b[2] = sigma_b2
    }
    
    # Validate P1, P2, P.DIFF & OR1, OR2, OR.DIFF
    parm1.arg.list = list(p1, p2, p.diff)
    parm1.args = unlist(lapply(parm1.arg.list, is.null))
    parm2.arg.list = list(or1, or2, or.diff)
    parm2.args = unlist(lapply(parm2.arg.list, is.null))
    if(sum(parm1.args) < 3 && sum(parm2.args) < 3){
      stop("Only one set of parameters may be supplied: Expected probabilities OR expected odds ratios")
    }
    if(sum(parm2.args) == 3 && sum(parm1.args) > 1){
      stop("At least two of the following terms must be specified: P1, P2, P.DIFF")
    }
    if(sum(parm1.args) == 3 && sum(parm2.args) > 1){
      stop("At least two of the following terms must be specified: OR1, OR2, OR.DIFF")
    }
    if(sum(parm1.args) == 0 && p.diff != abs(p1 - p2)){
      stop("At least one of the following terms has been misspecified: P1, P2, P.DIFF")
    }
    if(sum(parm2.args) == 0 && or.diff != abs(or1 - or2)){
      stop("At least one of the following terms has been misspecified: OR1, OR2, OR.DIFF")
    }
    
    # Validate ALPHA, METHOD, QUIET, ALL.SIM.DATA
    if(!is.numeric(alpha) || alpha < 0){
      stop("ALPHA", min0.warning)
    } else if(alpha > 1){
      stop("ALPHA must be a numeric value between 0 - 1")
    }
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
    
    # Calculate all expected probabilities/odds ratios (if they have not been specified)
    if(sum(parm2.args) == 3){
      if(is.null(p1)){
        p1 = abs(p.diff - p2)
      }
      if(is.null(p2)){
        p2 = abs(p1 - p.diff)
      }
      if(is.null(p.diff)){
        p.diff = abs(p1 - p2)
      }
    }
    if(sum(parm1.args) == 3){
      if(is.null(or1)){
        or1 = abs(or.diff - or2)
      }
      if(is.null(or2)){
        or2 = abs(or1 - or.diff)
      }
      if(is.null(or.diff)){
        or.diff = or1 - or2
      }
      p1 = or1 / (1 + or1)
      p2 = or2 / (1 + or2)
      p.diff = abs(p1 - p2)
    }
    
    # Calculate ICC1 (sigma_b / (sigma_b + pi^2/3))
    icc1 = mean(sapply(1:2, function(x) sigma_b[x] / (sigma_b[x] + pi^2 / 3)))
    
    # Create indicators for treatment group & cluster
    trt = c(rep(0, length.out = sum(nsubjects[1:nclusters[1]])), 
            rep(1, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
    clust = unlist(lapply(1:sum(nclusters), function(x) rep(x, length.out = nsubjects[x])))
    
    # Calculate log odds for each group
    logit.p1 = log(p1 / (1 - p1))
    logit.p2 = log(p2 / (1 - p2))
    
    # Set warnings to OFF
    # Note: Warnings will still be stored in 'warning.list'
    options(warn = -1)
    
    ### Create simulation loop
    while(sum(converge.vector == TRUE) != nsim){
      # Generate between-cluster effects for non-treatment and treatment
      randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b[1]))
      randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b[2]))
      
      # Create non-treatment y-value
      y0.intercept = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
      y0.linpred = y0.intercept + logit.p1
      y0.prob = expit(y0.linpred)
      y0 = unlist(lapply(y0.prob, function(x) stats::rbinom(1, 1, x)))
      if (length(table(y0))!=2){
        warning(print("y0 is completely seperated. Repeating the random draw 1 time."))
        randint.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b[1]))
        y0.intercept = unlist(lapply(1:nclusters[1], function(x) rep(randint.0[x], length.out = nsubjects[x])))
        y0.linpred = y0.intercept + logit.p1
        y0.prob = expit(y0.linpred)
        y0 = unlist(lapply(y0.prob, function(x) stats::rbinom(1, 1, x)))
      }

      # Create treatment y-value
      y1.intercept = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
      y1.linpred = y1.intercept + logit.p2
      y1.prob = expit(y1.linpred)
      y1 = unlist(lapply(y1.prob, function(x) stats::rbinom(1, 1, x)))
      if (length(table(y1))!=2){
        warning(print("y1 is completely seperated. Repeating the random draw 1 time."))
        randint.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b[2]))
        y1.intercept = unlist(lapply(1:nclusters[2], function(x) rep(randint.1[x], length.out = nsubjects[nclusters[1] + x])))
        y1.linpred = y1.intercept + logit.p2
        y1.prob = expit(y1.linpred)
        y1 = unlist(lapply(y1.prob, function(x) stats::rbinom(1, 1, x)))
      }

      
      
      # Create single response vector
      y = c(y0, y1)
      
      # Create and store data frame for simulated dataset
      sim.dat = data.frame(y = y, trt = trt, clust = clust)
      if(all.sim.data == TRUE){
        simulated.datasets = append(simulated.datasets, list(sim.dat))
      }
      # Calculate ICC2 ([P(Yij = 1, Yih = 1)] - pij * pih) / sqrt(pij(1 - pij) * pih(1 - pih))
      #icc2 = (mean(y0.prob) * mean(y1.prob) - p1*p2) / sqrt((p1 * (1 - p1)) * p2 * (1 - p2))
      icc2 = (mean(y0.prob) - p1) * (mean(y1.prob) - p2) / sqrt((p1 * (1 - p1)) * p2 * (1 - p2))
      # ^Equation above #11 (no number); Eldridge, Ukoumunne & Carlin, 2009 (p.386)
      icc2.vector = append(icc2.vector, icc2)
      
      # Calculate LMER.ICC (lmer: sigma_b / (sigma_b + sigma))
      if(irgtt==TRUE){
        my.mod <- lme4::lmer(y ~ trt + (0 + trt|clust), data = sim.dat, 
                             family = stats::binomial(link = 'logit'))
      } else {
      lmer.mod = lme4::lmer(y ~ trt + (1|clust), data = sim.dat)}
      lmer.vcov = as.data.frame(lme4::VarCorr(lmer.mod))[, 4]
      lmer.icc.vector = append(lmer.icc.vector, lmer.vcov[1] / (lmer.vcov[1] + lmer.vcov[2]))
      
      # Fit GLMM (lmer)
      if(method == 'glmm'){
        if(irgtt==TRUE){
          my.mod <- lme4::lmer(y ~ trt + (0 + trt|clust), data = sim.dat, 
                               family = stats::binomial(link = 'logit'))
        } else {
        my.mod = try(lme4::glmer(y ~ trt + (1|clust), data = sim.dat, 
                                 family = stats::binomial(link = 'logit')))
        }
        model.converge = try(my.mod)
        converge.ind = is.null(model.converge@optinfo$conv$lme4$messages)
        converge.vector = append(converge.vector, converge.ind)
        if(converge.ind == FALSE){
          model.id = paste0("Model ", length(converge.vector))
          warning.list[model.id] = list(model.converge@optinfo$conv$lme4$messages)
        }
        glmm.values = summary(my.mod)$coefficient
        est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
        se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
        stat.vector = append(stat.vector, glmm.values['trt', 'z value'])
        pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|z|)'])
      }
      # Fit GEE (geeglm)
      if(method == 'gee'){
        sim.dat = dplyr::arrange(sim.dat, clust)
        my.mod = geepack::geeglm(y ~ trt, data = sim.dat,
                                 family = stats::binomial(link = 'logit'), 
                                 id = clust, corstr = "exchangeable")
        gee.values = summary(my.mod)$coefficients
        est.vector = append(est.vector, gee.values['trt', 'Estimate'])
        se.vector = append(se.vector, gee.values['trt', 'Std.err'])
        stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
        pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
        converge.vector = append(converge.vector, TRUE)
      }
      
      # Update simulation progress information
      if(quiet == FALSE){
        # Print simulation start message
        if(length(est.vector) == 1){
          avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
          time.est = avg.iter.time * (nsim - 1) / 60
          hr.est = time.est %/% 60
          min.est = round(time.est %% 60, 0)
          message(paste0('Begin simulations :: Start Time: ', Sys.time(), 
                       ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
        }
        # Print simulation complete message
        if(sum(converge.vector == TRUE) == nsim){
          message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
        }
      }
      # Iterate progress bar
      prog.bar$update(sum(converge.vector == TRUE) / nsim)
      Sys.sleep(1/100)
      
      # Governor to prevent infinite non-convergence loop
      converge.ratio = sum(converge.vector == FALSE) / sum(converge.vector == TRUE)
      if(converge.ratio > 4.0 && converge.ratio != Inf){
        stop("WARNING! The number of non-convergent models exceeds the number of convergent models by a factor of 4. Consider reducing SIGMA_B")
      }
    }
    
    ## Output objects
    # Create object containing summary statement
    summary.message = paste0("Monte Carlo Power Estimation based on ", nsim, 
                             " Simulations: Simple Design, Binary Outcome\nNote: ", sum(converge.vector == FALSE), 
                             " additional models were fitted to account for non-convergent simulations.")
    
    # Create method object
    long.method = switch(method, glmm = 'Generalized Linear Mixed Model', 
                         gee = 'Generalized Estimating Equation')
    
    # Store model estimate output in data frame
    cps.model.est = data.frame(Estimate = as.vector(unlist(est.vector)),
                             Std.err = as.vector(unlist(se.vector)),
                             Test.statistic = as.vector(unlist(stat.vector)),
                             p.value = as.vector(unlist(pval.vector)), 
                             converge = as.vector(unlist(converge.vector)))
    cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
    
    # Calculate and store power estimate & confidence intervals
    # pval.data = subset(cps.model.est, converge == TRUE)
    pval.data = cps.model.est[cps.model.est$converge == TRUE,]
    pval.power = sum(pval.data[, 'sig.val']) / nrow(pval.data)
    power.parms = data.frame(power = round(pval.power, 3),
                             lower.95.ci = round(pval.power - abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                             upper.95.ci = round(pval.power + abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
    
    # Create object containing inputs
    p1.p2.or = round(p1 / (1 - p1) / (p2 / (1 - p2)), 3) 
    p2.p1.or = round(p2 / (1 - p2) / (p1 / (1 - p1)), 3) 
    inputs = t(data.frame('Non.Treatment' = c("probability" = p1, "odds.ratio" = p1.p2.or), 
                          'Treatment' = c("probability" = p2, 'odds.ratio' = p2.p1.or), 
                          'Difference' = c("probability" = p.diff, 'odds.ratio' = p2.p1.or - p1.p2.or)))

    # Create object containing group-specific cluster sizes
    cluster.sizes = list('Non.Treatment' = nsubjects[1:nclusters[1]], 
                         'Treatment' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
    
    # Create object containing number of clusters
    n.clusters = t(data.frame("Non.Treatment" = c("n.clust" = nclusters[1]), "Treatment" = c("n.clust" = nclusters[2])))
    
    # Create object containing estimated ICC values
    ICC = round(t(data.frame('P_h' = c('ICC' = icc1), 
                             'P_c' = c('ICC' = mean(icc2.vector)), 
                             'lmer' = c('ICC' = mean(lmer.icc.vector)))), 3)
    # Create object containing all ICC values
    # Note: P_h is a single calculated value. No vector to be appended.
    icc.list = data.frame('P_c' = icc2.vector, 
                          'lmer' = lmer.icc.vector)
    
    # Create object containing group-specific variance parameters
    var.parms = t(data.frame('Non.Treatment' = c('sigma_b' = sigma_b[1]), 
                             'Treatment' = c('sigma_b' = sigma_b[2])))
    
    # Check & governor for inclusion of simulated datasets
    # Note: If number of non-convergent models exceeds 5% of NSIM, override ALL.SIM.DATA and output all simulated data sets
    if(all.sim.data == FALSE && (sum(converge.vector == FALSE) < sum(converge.vector == TRUE) * 0.05)){
      simulated.datasets = NULL
    }
    
    # Create list containing all output (class 'crtpwr') and return
    complete.output = structure(list("overview" = summary.message, "nsim" = nsim, "power" = power.parms, "method" = long.method, "alpha" = alpha,
                                     "cluster.sizes" = cluster.sizes, "n.clusters" = n.clusters, "variance.parms" = var.parms, 
                                     "inputs" = inputs, "ICC" = ICC, "icc.list" = icc.list, "model.estimates" = cps.model.est, 
                                     "sim.data" = simulated.datasets, "warning.list" = warning.list), class = 'crtpwr')
    return(complete.output)
    }
