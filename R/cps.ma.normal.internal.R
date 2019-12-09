#' Runs and returns simulations and model fits for multi-arm 
#' cluster-randomized trials with continuous outcome.
#'
#' This function is called within \code{cps.ma.normal()} to generate the 
#' simulated data and return glmm or gee model fits according to the user's 
#' specifications. However, this function can also be called independently 
#' in order to return and examine the simulation model fits rather than the 
#' power summary returned in \code{cps.ma.normal()}. 
#' 
#' Users can modify a variety of parameters to suit the simulations to their
#' desired experimental situation. Users must specify the desired number 
#' of simulations, number of subjects per cluster, number of clusters 
#' per treatment arm, group means, within-cluster variance, and between-cluster 
#' variance; significance level, analytic method, progress updates, and 
#' simulated data set output may also be specified.
#'
#' As this function is not intended for the casual user, specifications to the 
#' arguments are somewhat less flexible than those supplied to 
#' \code{cps.ma.normal()}. Note that the str.nsubjects argument has a specific 
#' structure that is somewhat verbose but maximizes flexibility in the 
#' experimental design. 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param str.nsubjects Number of subjects per treatment group; accepts a list 
#' with one entry per arm. Each entry is a vector containing the number of 
#' subjects per cluster (required).
#' @param means Expected absolute treatment effect for each arm; accepts a 
#' vector of length \code{narms} (required).
#' @param sigma_sq Within-cluster variance; accepts a vector of length 
#' \code{narms} (required).
#' @param sigma_b_sq Between-cluster variance; accepts a vector of length 
#' \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects 
#' Model (GLMM) or Generalized Estimating Equation (GEE); accepts c('glmm', 
#' 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; 
#' default = FALSE.
#' @param seed Option to set.seed, default is NULL.
#' @param cores a string or numeric value indicating the number of cores to be 
#' used for parallel computing. When this option is set to NULL, no parallel 
#' computing is used.
#' @param poor.fit.override Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge.
#' @param low.power.override Option to override \code{stop()} if the power 
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the 
#' stop message. Default = FALSE. When TRUE, this check is ignored and the 
#' calculated power is returned regardless of value. 
#' @param tdist Logical; use t-distribution instead of normal distribution 
#' for simulation values, default = FALSE.
#' @param opt Option to fit with a different optimizer (using the package \code{optimx}).
#' @return A list with the following components:
#' \describe{
#'   \item{estimates}{List of \code{length(nsim)} containing gee- or glmm-fitted model 
#'   estimates.}
#'   \item{model comparisons}{Compares fitted model to a model for H0 using ML 
#'   (ANOVA).}
#'   \item{simulated datasets}{Produced when all.sim.data = TRUE, a list of 
#'   data frames, each containing three columns: "y" (Simulated response value),  
#'                   "trt" (Indicator for treatment group), and 
#'                   "clust" (Indicator for cluster).}
#'   \item{converged}{A logical vector of length \code{nsim} (when all.sim.data = TRUE) or a string (when all.sim.data = FALSE) 
#'   indicating the percent of model fits that produced a "singular fit" or "failed to converge" warning message;
#'           When a model fails to converge, failed.to.converge == FALSE, otherwise TRUE.}
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @examples 
#' \dontrun{
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' means.example <- c(22, 21, 21.5)
#' sigma_sq.example <- c(1, 1, 0.9)
#' sigma_b_sq.example <- c(0.1, 0.15, 0.1)
#' 
#' multi.cps.normal.models <- cps.ma.normal.internal (nsim = 100, 
#'                               str.nsubjects = nsubjects.example, 
#'                               means = means.example, 
#'                               sigma_sq = sigma_sq.example, 
#'                               sigma_b_sq = sigma_b_sq.example, 
#'                               alpha = 0.05,
#'                               quiet = FALSE, method = 'glmm', 
#'                               seed = 123, cores = "all",
#'                               low.power.override = FALSE,
#'                               poor.fit.override = FALSE, 
#'                               opt = "nlm")
#'                               }
#' @export

cps.ma.normal.internal <-  function(nsim = 1000, str.nsubjects = NULL,
                      means = NULL, sigma_sq = NULL, sigma_b_sq = NULL,
                      alpha = 0.05,
                      quiet = FALSE, method = 'glmm', 
                      all.sim.data = FALSE, 
                      seed = NA,
                      cores = "all",
                      poor.fit.override = FALSE,
                      low.power.override = FALSE,
                      tdist=FALSE,
                      opt = "optim"){

  # Create vectors to collect iteration-specific values
  simulated.datasets = list()
  
  # Create NCLUSTERS, NARMS, from str.nsubjects
  narms = length(str.nsubjects)
  nclusters = sapply(str.nsubjects, length)
  
  # Set start.time for progress iterator & initialize progress bar
  start.time = Sys.time()
  prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta", 
                                         total = nsim, clear = FALSE, width = 100)
  
  # This container keeps track of how many models failed to converge
  converge.vector <- rep(TRUE, nsim)
  
  # Create a container for the simulated.dataset and model output
  sim.dat = vector(mode = "list", length = nsim)
  model.values <- list()
  model.compare <- list()
  
  # option for reproducibility
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create indicators for treatment group & cluster for the sim.data output
  trt = list()
  for (arm in 1:length(str.nsubjects)) {
    trt[[arm]] = list()
    for (cluster in 1:length(str.nsubjects[[arm]])){
      trt[[arm]][[cluster]] = rep(arm, str.nsubjects[[arm]][[cluster]])
    }
  }
  clust = list()
  for(i in 1:sum(nclusters)){
    clust[[i]] <- lapply(seq(1, sum(nclusters))[i], 
                         function (x) {rep.int(x, unlist(str.nsubjects)[i])})
  }
  #Alert the user if using t-distribution
  if (tdist==TRUE){
    print("using t-distribution because tdist = TRUE")
  }
  #setup for parallel computing
  if (!is.null(cores)){
    ## Do computations with multiple processors:
    ## Number of cores:
    if (cores=="all"){nc <- parallel::detectCores()} else {nc <- cores}
    ## Create clusters
    cl <- parallel::makeCluster(rep("localhost", nc))
  }
  
  # Create simulation loop
  require(foreach)
  foreach::foreach(i=1:nsim, .packages = c("optimx")) %do% {
    sim.dat[[i]] = data.frame(y = NA, trt = as.factor(unlist(trt)), 
                              clust = as.factor(unlist(clust)))
    # Generate between-cluster effects for non-treatment and treatment
    if (tdist==TRUE){
      randint = mapply(function(n, df) stats::rt(n, df = df), 
                       n = nclusters, 
                       df = Inf)
    } else {
    randint = mapply(function(nc, s, mu) stats::rnorm(nc, mean = mu, sd = sqrt(s)), 
                                                      nc = nclusters, s = sigma_b_sq, 
                                                      mu = 0)
    }
    # Create y-value
    y.bclust <-  vector(mode = "numeric", length = length(unlist(str.nsubjects)))
    y.wclust <-  vector(mode = "list", length = narms)
    y.bclust <-  sapply(1:sum(nclusters), 
                      function(x) rep(unlist(randint)[x], length.out = unlist(str.nsubjects)[x]))
    for (j in 1:narms) {
      y.wclust[[j]] <-  lapply(str.nsubjects[[j]], function(x) stats::rnorm(x, mean = means[j], 
                                                                             sd = sqrt(sigma_sq[j])))
    }

    # Create data frame for simulated dataset
    y <- as.vector(unlist(y.bclust) + unlist(y.wclust))
    sim.dat[[i]][["y"]] <- y
  
    # Update simulation progress information
    if(quiet == FALSE){
      if(i == 1){
        message(paste0('Begin simulations :: Start Time: ', Sys.time()))
      }
      # Iterate progress bar
      prog.bar$update(i / nsim)
      Sys.sleep(1/100)
    }
    
    # trt and clust are re-coded as trt2 and clust2 to work nicely with lme. This can be changed later.
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      if (max(sigma_sq) != min(sigma_sq) & max(sigma_b_sq) != min(sigma_b_sq)){
        trt2 <- unlist(trt)
        clust2 <- unlist(clust)
      if (opt != "nlm" && opt != "nlminb" && opt != "auto"){
        stop("opt must be either nlm or nlminb for this model type.")
      }  
      my.mod <- try(nlme::lme(y~as.factor(trt2), random = ~1+as.factor(trt2)|clust2, 
                     weights = nlme::varIdent(form = ~1|as.factor(trt2)), 
                     method = "ML",
                     control = nlme::nlmeControl(opt = opt, niterEM = 100,
                                                 msMaxIter = 100)))
      model.values[[i]] <-  try(summary(my.mod)$tTable)
      # get the overall p-values (>Chisq)
      null.mod <- try(nlme::lme(y~1, random = ~1 + as.factor(trt2)|clust2, 
                            weights = nlme::varIdent(form = ~1|as.factor(trt2)), 
                                                   method="ML",
                            control=nlme::nlmeControl(opt = opt, niterEM = 100,
                                                      msMaxIter = 100)))
      converge.vector[i] <- ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE) 
      if (poor.fit.override==FALSE){
        if(sum(converge.vector == FALSE, na.rm = TRUE)>(nsim*.25)){stop("more than 25% of simulations are singular fit: check model specifications")
        }
      }
      }
      
      if (max(sigma_sq)==min(sigma_sq) & max(sigma_b_sq)!=min(sigma_b_sq)){
        my.mod <-  lmerTest::lmer(y~trt+(1+as.factor(trt)|clust), REML=FALSE,
                                  data = sim.dat[[i]])
        if (!isTRUE(class(my.mod) == "nlme")){
          if (i == 1){
            if (opt == "auto"){
              require("optimx")
              goodopt <- optimizerSearch(my.mod)

            } else {
              goodopt <- opt
            }
          }
        } 
        if (isTRUE(class(my.mod) == "nlme")){
          if (i == 1){
            if (opt == "auto"){
              goodopt <- "nlminb"
            } else {
              goodopt <- opt
            }
          }
        } 
        my.mod <-  lmerTest::lmer(y~trt+(1+as.factor(trt)|clust), REML=FALSE,
                                  data = sim.dat[[i]], 
                                  lme4::lmerControl(optimizer = goodopt))
        # get the overall p-values (>Chisq)
        null.mod <- update.formula(my.mod, y ~ 1 + (1+as.factor(trt)|clust))
        # option to stop the function early if fits are singular
        converge.vector[i] <- ifelse(isTRUE(any( grepl("singular", my.mod@optinfo$conv$lme4$messages))), FALSE, TRUE) 
        if (poor.fit.override==FALSE){
          if(sum(converge.vector == FALSE, na.rm = TRUE)>(nsim*.25)){stop("more than 25% of simulations are singular fit: check model specifications")}
        }
      }
      
      if (max(sigma_sq)!=min(sigma_sq) & max(sigma_b_sq)==min(sigma_b_sq)){
        trt2 <- unlist(trt)
        clust2 <- unlist(clust)
        if (opt != "nlm" && opt != "nlminb"){
          stop("opt must be either nlm or nlminb for this model type.")
        }  
        my.mod <- try(nlme::lme(y~as.factor(trt2), random = ~1|clust2, 
                            weights=nlme::varIdent(form = ~1|as.factor(trt2)), 
                            method="ML",
                            control=nlme::nlmeControl(opt = opt, niterEM = 100, 
                                                      msMaxIter = 100)))
        model.values[[i]] <-  try(summary(my.mod)$tTable)
        # get the overall p-values (>Chisq)
        null.mod <- try(nlme::lme(y~1, random=~1|clust2,  
                              weights=nlme::varIdent(form = ~1|as.factor(trt2)), 
                              method="ML",
                              control=nlme::nlmeControl(opt = opt, niterEM = 100, 
                                                        msMaxIter = 100)))
        converge.vector[i] <- ifelse(isTRUE(class(my.mod) == "try-error"), FALSE, TRUE) 
        if (poor.fit.override==FALSE){
          if(sum(converge.vector == FALSE, na.rm = TRUE)>(nsim*.25)){stop("more than 25% of simulations are singular fit: check model specifications")
          }
        }
      }
      
      if (max(sigma_sq)==min(sigma_sq) & max(sigma_b_sq)==min(sigma_b_sq)){
        my.mod <- lmerTest::lmer(y~trt+(1|clust), REML=FALSE, 
                                  data = sim.dat[[i]])
        if (i == 1){
          if (!isTRUE(class(my.mod) == "nlme")){
            if (opt == "auto"){
              require("optimx")
              goodopt <- optimizerSearch(my.mod)
            } else {
              goodopt <- opt
            }
          }
        }
        my.mod <-  lmerTest::lmer(y~trt+(1|clust), REML=FALSE, 
                                  data = sim.dat[[i]], 
                                  lme4::lmerControl(optimizer = goodopt))
        # get the overall p-values (>Chisq)
        null.mod <- update.formula(my.mod, y ~ 1 + (1|clust))
        # option to stop the function early if fits are singular
        converge.vector[i] <- ifelse(isTRUE(any(grepl("singular", my.mod@optinfo$conv$lme4$messages))), FALSE, TRUE) 
        if (poor.fit.override==FALSE){
          if(sum(converge.vector == FALSE, na.rm = TRUE)>(nsim*.25)){stop("more than 25% of simulations are singular fit: check model specifications")}
        }
      } 
      
      model.values[[i]] <-  summary(my.mod)
    }
  # Fit GEE (geeglm)
  if(method == 'gee'){
    data.holder = dplyr::arrange(sim.dat[[i]], clust)
    
    # Iterate progress bar
    prog.bar$update(i / nsim)
    Sys.sleep(1/100)
    
    my.mod = geepack::geeglm(y ~ trt, data = data.holder,
                             id = clust, corstr = "exchangeable")
    model.values[[i]] = summary(my.mod)
  }
    

  model.compare[[i]] <- try(anova(my.mod, null.mod))
  # stop the loop if power is <0.5
  if (low.power.override==FALSE){
    if (i > 50 & (i %% 10==0)){
    temp.power.checker <- matrix(unlist(model.compare[1:i]), ncol=6, nrow=i, 
                                 byrow=TRUE)
    sig.val.temp <-  ifelse(temp.power.checker[,6][1:i] < alpha, 1, 0)
    pval.power.temp <- sum(sig.val.temp)/i
    if (pval.power.temp < 0.5){
      stop(paste("Calculated power is < ", pval.power.temp, ", auto stop at simulation ", 
                 i, ". Set low.power.override==TRUE to run the simulations anyway.", sep = ""))
    }
    }
  }
    
    if(i == nsim){
      total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      hr.est = total.est %/% 3600
      min.est = total.est %/% 60
      sec.est = round(total.est %% 60, 0)
      message(paste0("Simulations Complete! Time Completed: ", Sys.time(), 
                     "\nTotal Runtime: ", hr.est, 'Hr:', min.est, 'Min:', sec.est, 'Sec'))
    }
  }

  ## Output objects
  if(all.sim.data == TRUE){
    complete.output.internal <-  list("estimates" = try(model.values),
                                      "model.comparisons" = try(model.compare),
                                      "converged" = converge.vector,
                                      "sim.data" = sim.dat)
  } else {
    complete.output.internal <-  list("estimates" = try(model.values),
                                      "model.comparisons" = try(model.compare),
                                      "converged" = converge.vector)
  }
  
  # turn off parallel computing
  if (!exists("cores", mode = "NULL")){
    parallel::stopCluster(cl)
  }
  
  return(complete.output.internal)
}