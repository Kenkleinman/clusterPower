#' Canned estimation functions for the power simulations.
#'
#' These functions are designed to be used by the power.sim.XXX functions as the
#' functions which estimate the treatment effect. They fit simple fixed and random
#' effects models and return the estimated treatment effect. These functions are
#' not designed to be called directly by the user.
#' 
#' @param dat observed data as a data.frame with columns named, "y", "trt" and
#'   "clust". "per" column is optional if period.var==0.
#' @param incl.period.effect indicator of whether to include a period effect
#' @param outcome.type one of "gaussian", "binomial", "poisson"
#' @param alpha the type I error rate
#' @return A numeric vector with the following three elements, in order: 
#'  [1] a point estimate for the treatment effect, [2] lower bound of (1-alpha)
#'  confidence interval, [3] lower bound of (1-alpha) confidence interval.
#'  
#' @export

fixed.effect <- function(dat, incl.period.effect, outcome.type, alpha) {
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

	if(incl.period.effect==0){
		fit <- glm(y ~ trt + clust,
			   data=dat,
			   family=outcome.type,
			   offset=offsets)
	} else {
		fit <- glm(y ~ trt + per + clust - 1,
			   data=dat,
			   family=outcome.type,
			   offset=offsets)
	}

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)$coef["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))

}

#' @rdname fixed.effect
#' @export

fixed.effect.cluster.level <- function(dat, incl.period.effect, outcome.type, alpha) {
	cols <- c("y", "at.risk.time")

	clust.dat <- aggregate(dat[,cols], FUN=sum,
			       list(clust=dat[,"clust"],
				    per=dat[,"per"],
				    trt=dat[,"trt"]))
	if(outcome.type=="poisson"){
		offsets <- log(clust.dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(clust.dat))
	}

	## set the formula
	if(incl.period.effect==0 & outcome.type!="binomial")
		form <- formula(y ~ trt + clust)
	if(incl.period.effect==0 & outcome.type=="binomial"){
		successes <- clust.dat[,"y"]
		failures <- clust.dat[,"at.risk.time"]-clust.dat[,"y"]
		form <- formula(cbind(successes, failures) ~ trt + clust)
	}
	if(incl.period.effect!=0 & outcome.type!="binomial")
		form <- formula(y ~ trt + per + clust - 1)
	if(incl.period.effect!=0 & outcome.type=="binomial") {
		successes <- clust.dat[,"y"]
		failures <- clust.dat[,"at.risk.time"]-clust.dat[,"y"]
		form <- formula(cbind(successes, failures) ~ trt + per + clust - 1)
	}
	
	fit <- glm(form, data=clust.dat, family=outcome.type, offset=offsets)

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)$coef["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))

}

#' @rdname fixed.effect
#' @export

random.effect <- function(dat, incl.period.effect, outcome.type, alpha) {
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

  if(incl.period.effect==0){
    if(outcome.type=="gaussian"){
      fit <- lme4::lmer(y ~ trt + (1|clust),
                  data=dat,
                  offset=offsets)
    }
    else {
      fit <- lme4::glmer(y ~ trt + (1|clust),
                   data=dat,
                   family=outcome.type,
                   offset=offsets)
    }} else {
      if(outcome.type=="gaussian"){
        fit <- lme4::lmer(y ~ trt + per + (1|clust) - 1,
                    data=dat,
                    offset=offsets)
      }
      else{
        fit <- lme4::glmer(y ~ trt + per + (1|clust) - 1,
                     data=dat,
                     family=outcome.type,
                     offset=offsets)
      }}

	n.clust <- length(unique(dat$clust))
	df <- n.clust - 2 ## based on k-2 in Donner & Klar p.118
	t <- qt(1 - alpha/2, df=df) * c(-1, 1)
	est <- coef(summary(fit))["trt", ]
	ci <- est["Estimate"] + t * est["Std. Error"]
	return(c(est["Estimate"], ci))
}

#' @rdname fixed.effect
#' @export

## based on Turner (1997) method C2
weighted.crossover.cluster.level <- function(dat, incl.period.effect, outcome.type, alpha) {
        if(outcome.type %in% c("poisson", "binomial"))
                stop("Currently, only Gaussian models are supported with this method.")
        if(unique(dat[,"per"])!=2)
                stop("This method only works for crosover studies with 2 periods.")

        ## takes the mean outcome by cluster-period
        clust.dat <- aggregate(dat[,"y"], FUN=mean,
                               list(clust=dat[,"clust"],
                                    per=as.numeric(dat[,"per"]),
                                    trt=dat[,"trt"]))
        ## makes a wide file with separate colums for treatment/control means
        clust.means <- reshape(clust.dat, idvar="clust", 
                               v.names=c("x", "per"), 
                               timevar="trt", direction="wide")
        clust.diffs <- clust.means[,"x.1"]-clust.means[,"x.0"]
        ws <- clust.means[,"per.1"]-clust.means[,"per.0"]
        
        clust.sizes <- table(dat$clust)
        wts <- clust.sizes/2 ## section 3.2.2 in Turner (1997)
        
        ## fit linear model
        fit <- lm(clust.diffs ~ ws, weights=wts)

        
}