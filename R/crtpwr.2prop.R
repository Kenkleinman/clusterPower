#' Power calculations for simple cluster randomized trials, binary outcome
#'
#' Compute the power of a simple cluster randomized trial with a binary outcome,
#' or determine parameters to obtain a target power.
#' 
#' This function implements the approach of Donner and Klar (2000). An estimate for the
#'  intracluster correlation coefficient (ICC) is used to calculate a design effect that 
#'  accounts for variance inflation due to clustering. There are several ways in
#'  which estimates for the ICC for a binary outcome can be calculated, as described by
#'  Wu, Crespi, and Wong (2012). 
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'   
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param nsubjects The mean of the cluster sizes.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param p1 The expected proportion in the treatment group.
#' @param p2 The proportion in the control group.
#' @param icc The intraclass correlation.
#' @param pooled Logical indicating if pooled standard error should be used.
#' @param p1inc Logical indicating if p1 is expected to be greater than p2.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 10 observations per cluster, no variation in cluster size, probability
#' # in condition 1 of .1 and condition 2 of .2, and icc = 0.1.
#' crtpwr.2prop(nsubjects=10 ,p1=.1, p2=.2, icc=.1)
#' # 
#' # The result, showimg nclusters of greater than 37, suggests 38 clusters per 
#' # condition should be used.
#' @references Donner A, Klar N. Design and Analysis of Cluster Randomization
#' Trials in Health Research. London; Arnold; 2000.
#' @references Wu S, Crespi CM, Wong WK. Comparison of Methods for Etimating Intraclass
#' Correlation Coefficient for Binary Responses in Cancer Prevention Cluster Randomized
#' Trials. Contemp Clin Trials. 2012; 33(5): 869-880. doi:10.1016/j.cct.2012.05.004 London: Arnold; 2000.
#' @export

crtpwr.2prop <- function(alpha = 0.05, power = 0.80,
                         nclusters = NA, nsubjects = NA, cv = 0,
                         p1 = NA, p2 = NA,
                         icc = NA, pooled = FALSE,
                         p1inc = TRUE,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, nclusters, nsubjects, cv, p1, p2, icc)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "cv", "p1", "p2", "icc")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'cv', 'p1', 'p2', or 'icc' must be NA.")
  }
  
  target <- neednames[needind]
  
  pwr <- quote({
    DEFF <- 1 + ((cv^2 + 1)*nsubjects - 1)*icc
    if (pooled) {
      p <- (p1+p2)/2
      sdd <- sqrt(p*(1 - p)*2*DEFF/(nclusters*nsubjects))
    } else {
      sdd <- sqrt((p1*(1-p1) + p2*(1-p2))*DEFF/(nclusters*nsubjects))
    }
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    pnorm(abs(p1 - p2)/sdd - zcrit, lower.tail = TRUE)# +
    #pnorm(-zcrit - (p1 - p2)/sdd, lower.tail = TRUE)
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
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol)$root
  }
  
  # calculate p1
  if (is.na(p1)) {
    if(p1inc){
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                    interval = c(p2 + 1e-7, 1 - 1e-7),
                    tol = tol)$root
    } else {
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                    interval = c(1e-7, p2 - 1e-7),
                    tol = tol)$root
    }
  }
  
  # calculate p2
  if (is.na(p2)) {
    if(p1inc){
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                    interval = c(1e-7, p1 - 1e-7),
                    tol = tol)$root
      
    } else {
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                    interval = c(p1 + 1e-7, 1 - 1e-7),
                    tol = tol)$root
    }
  }
  
  # calculate nsubjects
  if (is.na(nsubjects)) {
    nsubjects <- stats::uniroot(function(nsubjects) eval(pwr) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
  }
  
  # calculate cv
  if (is.na(cv)) {
    
    cv <- stats::uniroot(function(cv) eval(pwr) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
  }
  
  # calculate icc
  if (is.na(icc)){
    icc <- stats::uniroot(function(icc) eval(pwr) - power,
                   interval = c(1e-07, 1 - 1e-7),
                   tol = tol)$root
  }
  
  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample proportion power calculation: ", target, sep = "")
  # note <- "'nclusters' is the number of clusters in each group and 'nsubjects' is the number of individuals in each cluster."
  # structure(list(nclusters = nclusters, nsubjects = nsubjects, cv = cv,
  #                p1 = p1, p1dec = p1dec, p1inc = p1inc,
  #                p2 = p2, p2dec = p2dec, p2inc = p2inc,
  #                icc = icc, alpha = alpha, power = power,
  #                note = note, method = method),
  #           class = "power.htest")
  
}
