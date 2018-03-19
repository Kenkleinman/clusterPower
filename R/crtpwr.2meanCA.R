#' Power calculations for simple cluster randomized trials, continuous outcome with covariate adjustment
#'
#' Compute the power of a simple cluster randomized trial with a continuous outcome after adjusting for covariates,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{d}, \code{icc}, \code{varw}, \code{covdf}, and \code{rho_b} must be passed as \code{NA}.
#'   Note that \code{alpha} and \code{power} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'   
#' If \code{nsubjects} is a vector the values, \code{nclusters} will be recalculated
#'    using the values in \code{nsubjects}.
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param nsubjects The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param varw The within-cluster variation.
#' @param covdf The degrees of freedom used by the cluster-level covariates. If \code{covdf} is 0,
#'   the value of \code{rho_b} must be 0 as well.
#' @param rho_b The between-cluster residual correlation between the outcome and the covariates.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.8, 10 observations per cluster, a difference of 1 unit,  icc = 0.1,
#' # a variance of 5 units, adjusting for 1 covariate, and a between-cluster correlation of 0.7:
#' crtpwr.2meanCA(nsubjects=10 ,d=1, icc=.1, varw=5, covdf=1, rho_b=.7)
#' # 
#' # The result, showimg nclusters of greater than 12, suggests 13 clusters per condition should be used.
#' @references Moerbeek M, Teerenstra S. (2016) Power Analysis of Trials with Multilevel Data. 
#'   Boca Raton: CRC Press; 2016: 84-87
#' @export

crtpwr.2meanCA <- function(alpha = 0.05, power = 0.80, nclusters = NA,
                           nsubjects = NA, d = NA, icc = NA,
                           varw = NA, covdf = NA, rho_b = NA,
                           tol = .Machine$double.eps^0.25){
  
  # if nsubjects is a vector, 
  if(length(nsubjects) > 1){
    nvec <- nsubjects
    nsubjects <- mean(nvec)
    nsd <- stats::sd(nvec)
    cv <- nsd/nsubjects
    nclusters <- length(nvec)
  }
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  if(covdf == 0 && rho_b != 0){
    stop("Because 'covdf' is 0, 'rho_b' must be 0.")
  }
  
  # list of needed inputs
  needlist <- list(alpha, power, nclusters, nsubjects, d, icc, varw, covdf, rho_b)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "d", "icc", "varw", "covdf", "rho_b")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'd', 'icc', 'varw', 'covdf', and 'rho_b' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # design effect
    DEFF <- 1 + (nsubjects - 1)*icc*(1 - rho_b^2)
    
    # # variance inflation
    # # if nvec exists, calcuate exact relative efficiency
    # if (exists("nvec")) {
    #   if(method == "taylor"){
    #     a <- (1 - icc)/icc
    #     RE <- ((nsubjects + a)/nsubjects)*(sum((nvec/(nvec+a)))/nclusters) # exact relative efficiency
    #     VIF <- DEFF*RE
    #   } else{
    #     VIF <- 1 + ((cv^2 + 1)*nsubjects - 1)*icc
    #   }
    # } else if(!is.na(nsubjects)){
    #   if(method == "taylor"){
    #     L <- nsubjects*icc/DEFF
    #     REt <- 1/(1 - cv^2*L*(1 - L)) # taylor approximation
    #     VIF <- DEFF*REt
    #   } else {
    #     VIF <- 1 + ((cv^2 + 1)*nsubjects - 1)*icc
    #   }
    # }
    
    tcrit <- qt(alpha/2, 2*(nclusters - 1) - covdf, lower.tail = FALSE)
    
    ncp <- sqrt(nclusters*nsubjects/(2*DEFF)) * abs(d)/sqrt(varw)
    
    pt(tcrit, 2*(nclusters - 1) - covdf, ncp, lower.tail = FALSE) #+ pt(-tcrit, 2*(nclusters - 1), ncp, lower.tail = TRUE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(pwr) - power,
                            interval = c(1e-10, 1 - 1e-10),
                            tol = tol)$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(pwr)
  }
  
  # calculate nclusters
  if (is.na(nclusters)) {
    nclusters <- stats::uniroot(function(nclusters) eval(pwr) - power,
                        interval = c(2 + covdf + 1e-10, 1e+7),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate nsubjects
  if (is.na(nsubjects)) {
    nsubjects <- stats::uniroot(function(nsubjects) eval(pwr) - power,
                        interval = c(2 + 1e-10, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate d
  if (is.na(d)) {
    d <- stats::uniroot(function(d) eval(pwr) - power,
                        interval = c(1e-10, 1e+07),
                        tol = tol)$root
  }
  
  # calculate icc
  if (is.na(icc)){
    icc <- stats::uniroot(function(icc) eval(pwr) - power,
                          interval = c(1e-07, 1 - 1e-07),
                          tol = tol)$root
  }
  
  # calculate varw
  if (is.na(varw)) {
    varw <- stats::uniroot(function(varw) eval(pwr) - power,
                           interval = c(1e-10, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  # calculate covdf
  if (is.na(covdf)) {
    covdf <- stats::uniroot(function(covdf) eval(pwr) - power,
                            interval = c(1e-10, 1e+07),
                            tol = tol, extendInt = "downX")$root
  }
  
  # calculate rho_b
  if (is.na(rho_b)) {
    rho_b <- stats::uniroot(function(rho_b) eval(pwr) - power,
                            interval = c(1e-07, 1 - 1e-07),
                            tol = tol)$root
  }
  
  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample t-test power calculation: ", target, sep = "")
  # note <- "'nclusters' is the number of clusters in each group and 'nsubjects' is the number of individuals in each cluster."
  # structure(list(alpha = alpha, power = power, nclusters = nclusters, nsubjects = nsubjects, cv = cv, d = d,
  #                icc = icc, varw = varw, note = note, method = method),
  #           class = "power.htest")
  
}
