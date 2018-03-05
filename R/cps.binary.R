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
#' ; significance level, analytic method, and whether or not progress updates should 
#' be displayed while the function is running.
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param n Number of subjects per cluster; accepts integer (required). 
#' @param m Number of clusters per treatment group; accepts integer (required).
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
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param alpha Significance level. Default = 0.05
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default is FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{sim.data}{Data frame with columns "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", "is.signif" (Is p-value less than alpha?)}
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
#' }
#' 
#' @author Alexander R. Bogdan
#' 
#' @examples 
#' \dontrun{
#' my.binary.sim = cps.binary(nsim = 100, n = 50, m = 6, p1 = 0.4, p2 = 0.2, sigma_b = 100,
#'                     alpha = 0.05, method = 'glmm', quiet = FALSE)
#' }
#'
#' @export

# Define function
cps.binary = function(nsim = NULL,m = NULL, n = NULL, p.diff = NULL,
                        p1 = NULL, p2 = NULL, or1 = NULL, or2 = NULL, or.diff = NULL, 
                        sigma_b = NULL, sigma_b2 = NULL, 
                        alpha = 0.05, method = 'glmm', quiet = FALSE){
    # Create vectors to collect iteration-specific values
    est.vector = NULL
    se.vector = NULL
    stat.vector = NULL
    pval.vector = NULL
    start.time = Sys.time()
    
    # Define wholenumber function
    is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    # Define expit function
    expit = function(x)  1 / (1 + exp(-x))
    
    # Validate NSIM, N, M, SIGMA_B, ALPHA
    sim.data.arg.list = list(nsim, n, m, sigma_b)
    sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
    if(sum(sim.data.args) > 0){
      stop("NSIM, N, M & SIGMA_B must all be specified. Please review your input values.")
    }
    min1.warning = " must be an integer greater than or equal to 1"
    if(!is.wholenumber(nsim) || nsim < 1){
      stop(paste0("NSIM", min1.warning))
    }
    if(!is.wholenumber(n) || n < 1){
      stop(paste0("N", min1.warning))
    }
    if(!is.wholenumber(m) || m < 1){
      stop(paste0("M", min1.warning))
    }
    if(length(n) > 2){
      stop("N can only be a vector of length 1 (equal # of clusters per group) or 2 (unequal # of clusters per group)")
    }
    
    min0.warning = " must be a numeric value greater than 0"
    if(!is.numeric(sigma_b) || sigma_b <= 0){
      stop("SIGMA_B", min0.warning)
    }
    if(!is.null(sigma_b2) && sigma_b2 <= 0){
      stop("SIGMA_B2", min0.warning)
    }
    
    if(!is.numeric(alpha) || alpha < 0){
      stop("ALPHA", min0.warning)
    } else if(alpha > 1){
      stop("ALPHA must be a numeric value between 0 - 1")
    }
    
    # Set cluster sizes for treatment arm (if not already specified)
    if(length(n) == 1){
      n[2] = n[1]
    }
    
    # Set sample sizes for each cluster (if not already specified)
    if(length(m) == 1){
      m[1:sum(n)] = m
    } 
    if(n[1] == n[2] && length(m) == n[1]){
      m = rep(m, 2)
    }
    if(length(n) == 2 && length(m) != 1 && length(m) != sum(n)){
      stop("A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for M")
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
      stop("At least one of the following terms has be misspecified: P1, P2, P.DIFF")
    }
    if(sum(parm2.args) == 0 && or.diff != abs(or1 - or2)){
      stop("At least one of the following terms has be misspecified: OR1, OR2, OR.DIFF")
    }
    
    # Simulation parameters
    if(sum(parm2.args) == 3){
      if(is.null(p1)){
        p1 = abs(p.diff - p2)
      }
      if(is.null(p2)){
        p2 = abs(p1 - p.diff)
      }
    }
    if(sum(parm1.args) == 3){
      if(is.null(or1)){
        or1 = abs(or.diff - or2)
      }
      if(is.null(or2)){
        or2 = abs(or1 - or.diff)
      }
      p1 = or1 / (1 + or1)
      p2 = or2 / (1 + or2)
    }
    
    # Validate METHOD, QUIET
    if(!is.element(method, c('glmm', 'gee'))){
      stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model) 
           or 'gee'(Generalized Estimating Equation)")
    }
    if(!is.logical(quiet)){
      stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
    }
    
    # Set between-cluster variances
    if(is.null(sigma_b2)){
      sigma_b[2] = sigma_b
    }else{
      sigma_b[2] = sigma_b2
    }
    
    # Create simulation loop
    for(i in 1:nsim){
      # Generate simulated data
      # Create indicators for treatment group & cluster
      trt = c(rep(0, length.out = sum(m[1:n[1]])), rep(1, length.out = sum(m[(n[1]+1):(n[1]+n[2])])))
      clust = unlist(lapply(1:sum(n), function(x) rep(x, length.out = m[x])))
      
      # Generate between-cluster effects for non-treatment and treatment
      randint.0 = stats::rnorm(n[1], mean = 0, sd = sqrt(sigma_b[1]))
      randint.1 = stats::rnorm(n[2], mean = 0, sd = sqrt(sigma_b[2]))
      
      # Create beta
      beta = log((p1 / (1 - p1)) / (p2 / (1 - p2)))
      
      # Create non-treatment y-value
      y0.intercept = unlist(lapply(1:n[1], function(x) rep(randint.0[x], length.out = m[x])))
      y0.linpred = y0.intercept + log(p1 / (1 - p1))
      y0.prob = expit(y0.linpred)
      y0 = unlist(lapply(y0.prob, function(x) rbinom(1, 1, x)))

      # Create treatment y-value
      y1.intercept = unlist(lapply(1:n[2], function(x) rep(randint.1[x], length.out = m[n[1]+x])))
      y1.linpred = y1.intercept + log(p2 / (1 - p2)) + beta
      y1.prob = expit(y1.linpred)
      y1 = unlist(lapply(y1.prob, function(x) rbinom(1, 1, x)))
      
      # Create single response vector
      y = c(y0,y1)
      
      # Create data frame for simulated dataset
      sim.dat = data.frame(y.resp = y, trt = trt, clust = clust)
      
      # Fit GLMM (lmer)
      if(method == 'glmm'){
        my.mod = lme4::glmer(y.resp ~ trt + (1|clust), data = sim.dat, family = binomial(link = 'logit'))
        #suppressWarnings(lme4::glmer(y.resp ~ trt + (1|clust), data = sim.dat, family = binomial(link = 'logit')))
        glmm.values = summary(my.mod)$coefficient
        est.vector = append(est.vector, glmm.values['trt', 'Estimate'])
        se.vector = append(se.vector, glmm.values['trt', 'Std. Error'])
        stat.vector = append(stat.vector, glmm.values['trt', 'z value'])
        pval.vector = append(pval.vector, glmm.values['trt', 'Pr(>|z|)'])
      }
      # Fit GEE (geeglm)
      if(method == 'gee'){
        sim.dat = dplyr::arrange(sim.dat, clust)
        my.mod = geepack::geeglm(y.resp ~ trt, data = sim.dat,
                                 family = binomial(link = 'logit'), 
                                 id = clust, corstr = "exchangeable")
        gee.values = summary(my.mod)$coefficients
        est.vector = append(est.vector, gee.values['trt', 'Estimate'])
        se.vector = append(se.vector, gee.values['trt', 'Std.err'])
        stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
        pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
      }
      
      if(quiet == FALSE){
        if(i == 1){
          avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
          time.est = avg.iter.time * (nsim - 1) / 60
          hr.est = time.est %/% 60
          min.est = round(time.est %% 60, 0)
          message(paste0('Begin simulations :: Start Time: ', Sys.time(), ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
        }
        else if(i == nsim){
          message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
        } 
        else if(i %% 10 == 0){
          time.est = avg.iter.time * (nsim - i) / 60 
          hr.est = time.est %/% 60
          min.est = round(time.est %% 60, 0)
          min.est = ifelse(min.est == 0, '<1', min.est)
          message(paste0('Progress: ', i / nsim * 100, '% complete :: Estimated time remaining: ', hr.est, 'Hr:', min.est, 'Min'))
        }
      }
    }
    cps.sim.dat = data.frame(estimates = as.vector(unlist(est.vector)), 
                             stderrs = as.vector(unlist(se.vector)),
                             test.stat = as.vector(unlist(stat.vector)),
                             pvals = as.vector(unlist(pval.vector)))
    cps.sim.dat[, 'sig.vals'] = ifelse(cps.sim.dat[, 'pvals'] < alpha, 1, 0)
    pval.power = sum(cps.sim.dat[, 'sig.vals']) / length(cps.sim.dat[, 'sig.vals'])
    power.parms = data.frame(power = round(pval.power, 3), 
                             lower.95.ci = round(pval.power - abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3), 
                             upper.95.ci = round(pval.power + abs(qnorm(alpha/2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
    complete.output = list("sim.data" = cps.sim.dat, "power" = power.parms)
    return(complete.output)
    }
