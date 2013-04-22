## estimation functions should have
## INPUTS
##   ~ dat = data as a data.frame
##   ~ incl.period.effect = indicator of whether to include a period effect
##   ~ outcome.type = one of "gaussian", "binomial", "poisson"
##   ~ alpha = the type 1 error rate
## OUTPUTS
##   ~ a vector with three elements, in order:
##       [1] a point estimate for the treatment effect
##       [2] lower bound of (1-alpha) confidence interval
##       [3] lower bound of (1-alpha) confidence interval


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


random.effect <- function(dat, incl.period.effect, outcome.type, alpha) {
	require(lme4)
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

	if(incl.period.effect==0){
		fit <- glmer(y ~ trt + (1|clust),
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	} else {
		fit <- glmer(y ~ trt + per + (1|clust) - 1,
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	}

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)@coefs["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))
}
