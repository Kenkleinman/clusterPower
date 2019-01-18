# Validation functions
validateVariance <- function(x){
  warning("FIXME: not actually validating variance yet")
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


##FIXME: TO DO
# 11. testthat tests
# 2. update the example/ man text
# 3. input validation
# 4. make validateVariance fxn in "validation" file
# 5. make validate nsubjects fxn in validation file
# 6. make a wrapper function that takes ICC or sigma/sigma_b, nclusters?, narms?, formats output
# 7. make a seperate function for taking ICC, sigma, sigma_b
# 9. write some usage examples
# 13. Make sure man page for the wrapper also notes that the responsibility 
# for correcting for multiple testing lies with the user.
# 14. Must be able to set the seed on the simulation methods.
# 15. set.seed() option in the wrapper


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
#' 
#' @author Alexandria C. Sakrejda
#' @author Alexander R. Bogdan
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per treatment group; accepts a list with one entry per arm. 
#' Each entry is a vector containing the number of subjects per cluster (required).
#' @param means Expected absolute treatment effect for each arm; accepts a vector of length \code{narms} (required).
#' @param sigma Within-cluster variance; accepts a vector of length \code{narms} (required).
#' @param sigma_b Between-cluster variance; accepts a vector of length \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{model.values}{List of length(nsim) containing gee- or glmm-fitted the model summaries.
#'   Note: the responsibility for correcting for multiple testing lies with the user.}
#'   \item{sim.data}{List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)}
#' }
#' 
#' 
#' 
#' @examples 
#' \dontrun{
#' 
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' means.example <- c(30, 21, 53)
#' sigma.example <- c(1, 1, 0.9)
#' sigma_b.example <- c(0.1, 0.15, 0.1)
#' 
#' normal.ma.rct <- cps.ma.normal.internal(nsim = 100, nsubjects = nsubjects.example, 
#'                                        means = means.example, sigma = sigma.example, 
#'                                        sigma_b = sigma_b.example, alpha = 0.05, 
#'                                        quiet = FALSE, method = 'glmm', 
#'                                        all.sim.data = FALSE, seed = 123)
#' }
#' 
#' @export

cps.ma.normal.internal = function(nsim = NULL, nsubjects = NULL,
                      means = NULL, sigma = NULL, sigma_b = NULL,
                      alpha = 0.05,
                      quiet = FALSE, method = 'glmm', 
                      all.sim.data = FALSE, 
                      seed=NULL){

  # Create vectors to collect iteration-specific values
  simulated.datasets = list()
  
  # Create NCLUSTERS, NARMS, from NSUBJECTS
  narms = length(nsubjects)
  nclusters = sapply(nsubjects, length)
  
  # validation goes here
  
  # Set start.time for progress iterator & initialize progress bar
  start.time = Sys.time()
  prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta", 
                                         total = nsim, clear = FALSE, width = 100)
  prog.bar$tick(0)
  
  # Create indicators for treatment group & cluster
  trt = list()
  for (arm in 1:length(nsubjects)){
    trt[[arm]] = list()
    for (cluster in 1:length(nsubjects[[arm]])){
      trt[[arm]][[cluster]] = rep(arm, nsubjects[[arm]][[cluster]])
    }
  }
  clust = list()
  for (arm in 1:length(nsubjects)){
  clust[[arm]] = list()
    for (cluster in 1:length(nsubjects[[arm]])){
    clust[[arm]][[cluster]] = rep(cluster, nsubjects[[arm]][[cluster]])
    }
  }
  
  # This container keeps track of how many models failed to converge
  fail <- rep(NA, nsim)
  
  # Create a container for the simulated.dataset and glmm output
  sim.dat = vector(mode = "list", length = nsim)
  model.values <- list()
  
  # option for reproducibility
  set.seed(seed=seed)
  
  # Create simulation loop
  for(i in 1:nsim){
    sim.dat[[i]] = data.frame(y = NA, trt = as.factor(unlist(trt)), clust = as.factor(unlist(clust)))
    # Generate between-cluster effects for non-treatment and treatment
    randint = mapply(function(nc, s, mu) stats::rnorm(nc, mean = mu, sd = sqrt(s)), 
                                                      nc = nclusters, s = sigma_b, 
                                                      mu = means)
    # Create y-value
    y.bclust = vector(mode = "list", length = narms)
    y.wclust = vector(mode = "list", length = narms)
    y = vector(mode = "list", length = narms)
    for (j in 1:narms){
      y.bclust[[j]] = sapply(1:nclusters[j], function(x) rep(randint[[j]][x], length.out = nsubjects[[j]][x]))
      y.wclust[[j]] = lapply(nsubjects[[j]], function(x) stats:: rnorm(x, mean = randint[[j]], sd = sqrt(sigma[j])))
      y[[j]] = unlist(y.bclust[[j]]) + unlist(y.wclust[[j]])
    }
    # Create data frame for simulated dataset
    sim.dat[[i]][["y"]] = unlist(y)
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + (1|clust), data = sim.dat[[i]])
      model.values[[i]] = summary(my.mod)
      simulated.datasets[[i]] = sim.dat[[i]]
      }
  
  # Fit GEE (geeglm)
  if(method == 'gee'){
    data.holder = dplyr::arrange(sim.dat[[i]], clust)
    my.mod = geepack::geeglm(y ~ trt, data = data.holder,
                             id = clust, corstr = "exchangeable")
    model.values[[i]] = summary(my.mod)
 }
    # stop the function early if fits are singular
  fail[i] <- ifelse(any( grepl("singular", my.mod@optinfo$conv$lme4$messages) )==TRUE, 1, 0) 
  if(sum(fail, na.rm = TRUE)>=(nsim*.25)){stop("more than 25% of simulations are singular fit: check model specifications")}
  
  # Update simulation progress information
  if(quiet == FALSE){
    if(i == 1){
      avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      time.est = avg.iter.time * (nsim - 1) / 60
      hr.est = time.est %/% 60
      min.est = round(time.est %% 60, 0)
      message(paste0('Begin simulations :: Start Time: ', Sys.time(), 
                     ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
    }
    
    # Iterate progress bar
    prog.bar$update(i / nsim)
    Sys.sleep(1/100)
    
    if(i == nsim){
      total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      hr.est = total.est %/% 3600
      min.est = total.est %/% 60
      sec.est = round(total.est %% 60, 0)
      message(paste0("Simulations Complete! Time Completed: ", Sys.time(), 
                     "\nTotal Runtime: ", hr.est, 'Hr:', min.est, 'Min:', sec.est, 'Sec'))
    }
  }
  } 
  
  ## Output objects
  if(all.sim.data == TRUE){
    complete.output = list("estimates" = model.values,
                          "sim.data" = simulated.datasets,
                          "failed.to.converge"= fail)
  } else {
    complete.output = model.values
  }
  return(complete.output)
}