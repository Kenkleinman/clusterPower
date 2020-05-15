#' Model fits for simulations for multi-arm designs with dichotomous outcome.
#'
#' Generally called from \code{cps.ma.binary()}, this function uses iterative 
#' simulations to model significance of treatment effects for cluster-randomized 
#' controlled trials. Users can modify a variety of parameters to suit the 
#' simulations to their desired experimental situation. 
#' 
#' This function can be called directly in order to give the user access to the simulated 
#' model fits in addition to the simulated data, the latter of which can also be accessed 
#' here or using the function \code{cps.ma.binary()}. For the power estimates, use 
#' \code{cps.ma.binary()}.
#' 
#' Users (or the wrapper function) must specify the desired number of 
#' simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, group proportions, and 
#' between-cluster variance; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param str.nsubjects Number of subjects per treatment group; accepts a list with one entry per arm. 
#' Each entry is a vector containing the number of subjects per cluster (required).
#' @param probs Expected probability of outcome for each arm; accepts a vector of length \code{narms} (required).
#' @param sigma_b_sq Between-cluster variance; accepts a vector of length \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects 
#' Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 
#' 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; 
#' default = FALSE.
#' @param seed Option to set.seed. Default is NULL.
#' @param poor.fit.override Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge.
#' @param low.power.override Option to override \code{stop()} if the power 
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the 
#' stop message. Default = FALSE. When TRUE, this check is ignored and the 
#' calculated power is returned regardless of value. 
#' @param tdist Logical; use t-distribution instead of normal distribution for 
#' simulation values, default = FALSE.
#' @param cores A string ("all") NA, or numeric value indicating the number of cores to be used for parallel computing. 
#' When this option is set to NA, no parallel computing is used.
#' @param opt Option to fit with a different optimizer (using the package \code{optimx}). Default is 'optim'.
#' @param optmethod User-specified optimizer methods available for the optimizer specified in \code{opt} option.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item List of length(nsim) containing gee- or glmm-fitted the model summaries.
#'   \item Compares fitted model to a model for H0 using ML (anova).
#'   \item List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)
#'   \item A vector of logical values \code{nsim}. 
#'           When a model fails to converge, FALSE, otherwise TRUE.
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' nsubjects.example <- list(c(200,200,200,250), c(150, 200, 200, 210), c(170, 200, 210))
#' probs.example <- c(0.25, 0.15, 0.3)
#' sigma_b_sq.example <- c(0.01, 0.01, 0.01)
#' 
#' bin.ma.rct <- cps.ma.binary.internal (nsim = 10, 
#'                             str.nsubjects = nsubjects.example, 
#'                             probs = probs.example,
#'                             sigma_b_sq = sigma_b_sq.example, 
#'                             alpha = 0.05, all.sim.data = FALSE,
#'                             poor.fit.override = TRUE, 
#'                             seed = 123, cores="all",
#'                             opt = "bobyqa",
#'                             optmethod = "Nelder-Mead") 
#' }
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @export
cps.ma.binary.internal <-
  function(nsim = 1000,
           str.nsubjects = NULL,
           probs = NULL,
           sigma_b_sq = NULL,
           alpha = 0.05,
           quiet = FALSE,
           method = 'glmm',
           all.sim.data = FALSE,
           seed = NA,
           poor.fit.override = FALSE,
           low.power.override = FALSE,
           tdist = FALSE,
           cores = cores,
           opt = opt,
           optmethod = optmethod) {
    # Create vectors to collect iteration-specific values
    simulated.datasets = list()
    
    # Create NCLUSTERS, NARMS, from str.nsubjects
    narms = length(str.nsubjects)
    nclusters = sapply(str.nsubjects, length)
    
    # This container keeps track of how many models converged
    converged <- rep(NA, nsim)
    
    # Create a container for the simulated.dataset and model output
    sim.dat = vector(mode = "list", length = nsim)
    model.values <- list()
    model.compare <- list()
    
    # option for reproducibility
    if (!is.na(seed)) {
      set.seed(seed = seed)
    }
    
    # Create indicators for treatment group & cluster for the sim.data output
    trt1 = list()
    for (arm in 1:length(str.nsubjects)) {
      trt1[[arm]] = list()
      for (cluster in 1:length(str.nsubjects[[arm]])) {
        trt1[[arm]][[cluster]] = rep(arm, str.nsubjects[[arm]][[cluster]])
      }
    }
    clust1 = list()
    for (i in 1:sum(nclusters)) {
      clust1[[i]] <- lapply(seq(1, sum(nclusters))[i],
                            function (x) {
                              rep.int(x, unlist(str.nsubjects)[i])
                            })
    }
    
    # Calculate log odds for each group
    logit.p <- list()
    for (i in 1:length(probs)) {
      logit.p[[i]] = log(probs[i] / (1 - probs[i]))
    }
    logit.p <- unlist(logit.p)
    
    #Alert the user if using t-distribution
    if (tdist == TRUE) {
      print("using t-distribution because tdist = TRUE")
    }
    #make the simulated data
    trt <-  as.factor(unlist(trt1))
    clust <- as.factor(unlist(clust1))
    if (length(trt) != length(clust)) {
      stop("trt and clust are not the same length, see line 134")
    }
    sim.dat <- matrix(nrow = length(clust), ncol = nsim)
    # function to produce the simulated data
    make.sim.dat <- function(tdist = tdist,
                             logit.p = logit.p,
                             nclusters = nclusters,
                             sigma_b_sq = sigma_b_sq,
                             str.nsubjects = str.nsubjects) {
      # Generate between-cluster effects for non-treatment and treatment
      if (tdist == TRUE) {
        randint = mapply(function(n, df)
          stats::rt(n, df = df),
          n = nclusters,
          df = Inf)
      } else {
        randint = mapply(
          function(nc, s, mu)
            stats::rnorm(nc, mean = mu, sd = sqrt(s)),
          nc = nclusters,
          s = sigma_b_sq,
          mu = 0
        )
      }
      if (typeof(randint) == "list") {
        randint.holder <- list()
        for (j in 1:length(logit.p)) {
          randint.holder[[j]] <- logit.p[j] + randint[[j]]
        }
        randintrandint <- sapply(randint.holder, clusterPower::expit)
      } else {
        randint.holder <-
          matrix(nrow = nclusters[1], ncol = length(logit.p))
        for (j in 1:length(logit.p)) {
          randint.holder[, j] <- logit.p[j] + randint[, j]
        }
        randintrandint <- clusterPower::expit(randint.holder)
      }
      # Create y-value
      y.intercept <-  vector(mode = "numeric",
                             length = length(unlist(str.nsubjects)))
      y.intercept <-  sapply(1:sum(nclusters),
                             function(x)
                               rep(unlist(randintrandint)[x],
                                   length.out = unlist(str.nsubjects)[x]))
      y <-
        sapply(unlist(y.intercept), function(x)
          stats::rbinom(1, 1, x))
      return(y)
    }
    sim.dat <-
      replicate(
        nsim,
        make.sim.dat(
          tdist = tdist,
          logit.p = logit.p,
          nclusters = nclusters,
          sigma_b_sq = sigma_b_sq,
          str.nsubjects = str.nsubjects
        )
      )
    
    require(foreach)
    `%fun%` <- `%dopar%`
    if (is.na(cores)) {
      `%fun%` <- `%do%`
    }
    
    #setup for parallel computing
    ## Do computations with multiple processors:
    ## Number of cores:
    if (!is.na(cores)) {
      if (cores == "all") {
        nc <- parallel::detectCores()
      } else {
        nc <- cores
      }
      ## Create clusters and initialize the progress bar
      cl <- parallel::makeCluster(rep("localhost", nc), type = "SOCK")
      doSNOW::registerDoSNOW(cl)
    }
    pb <- txtProgressBar(min = 1, max = nsim, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    if (method == "glmm") {
      # Update simulation progress information
      if (quiet == FALSE) {
        sim.start <- Sys.time()
        lme4::glmer(sim.dat[, 1] ~ trt + (1 |
                                            clust),
                    family = stats::binomial(link = 'logit'))
        avg.iter.time = as.numeric(difftime(Sys.time(), sim.start, units = 'secs'))
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
        # initialize progress bar
        if (is.na(cores)) {
          prog.bar =  progress::progress_bar$new(
            format = "(:spin) [:bar] :percent eta :eta",
            total = 5,
            clear = FALSE,
            show_after = 0
          )
          prog.bar$tick(0)
        }
      }
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(1 / 5)
        Sys.sleep(1 / 100)
      }
      # Create simulation loop
      require(doParallel)
      require(foreach)
      if (!is.na(cores) & quiet == FALSE) {
        message("Fitting models")
      }
      
      my.mod <- foreach::foreach(
        i = 1:nsim,
        .options.snow = opts,
        .packages = c("lme4", "optimx"),
        .inorder = FALSE
      ) %fun% {
        lme4::glmer(
          sim.dat[, i] ~ trt + (1 | clust),
          family = stats::binomial(link = 'logit'),
          control = lme4::glmerControl(
            optimizer = opt,
            calc.derivs = TRUE,
            optCtrl = list(
              maxfun = 2e4,
              method = optmethod,
              starttests = TRUE,
              kkt = FALSE
            )
          )
        )
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(3 / 5)
        Sys.sleep(1 / 100)
      }
      
      # option to stop the function early if fits are singular
      for (i in 1:nsim) {
        converged[i] <-
          ifelse(is.null(my.mod[[i]]@optinfo$conv$lme4$messages),
                 TRUE,
                 FALSE)
      }
      
      if (poor.fit.override == FALSE) {
        if (sum(unlist(converged), na.rm = TRUE) < (nsim * .75)) {
          stop("more than 25% of simulations are singular fit: check model specifications")
        }
      }
      
      # stop the loop if power is <0.5
      if (low.power.override == FALSE) {
        if (i > 50 & (i %% 10 == 0)) {
          temp.power.checker <-
            matrix(
              unlist(model.compare[1:i]),
              ncol = 3,
              nrow = i,
              byrow = TRUE
            )
          sig.val.temp <-
            ifelse(temp.power.checker[, 3][1:i] < alpha, 1, 0)
          pval.power.temp <- sum(sig.val.temp) / i
          if (pval.power.temp < 0.5) {
            stop(
              paste(
                "Calculated power is < ",
                pval.power.temp,
                ", auto stop at simulation ",
                i,
                ". Set low.power.override==TRUE to run the simulations anyway.",
                sep = ""
              )
            )
          }
        }
      }
      
      if (!is.na(cores) & quiet == FALSE) {
        message("\r Performing null model comparisons")
      }
      # get the overall p-values (>Chisq)
      model.compare <- foreach::foreach(
        i = 1:nsim,
        .options.snow = opts,
        .packages = c("car", "optimx"),
        .inorder = FALSE
      ) %fun% {
        car::Anova(my.mod[[i]], type = "II")
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(4 / 5)
        Sys.sleep(1 / 100)
      }
      
      # get the model summaries
      if (!is.na(cores) & quiet == FALSE) {
        message("\r Retrieving model summaries")
      }
      model.values <-
        foreach::foreach(
          i = 1:nsim,
          .options.snow = opts,
          .packages = "car",
          .inorder = FALSE
        ) %fun% {
          summary(my.mod[[i]])
        }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(5 / 5)
        Sys.sleep(1 / 100)
      }
      # turn off parallel computing
      if (!is.na(cores)) {
        #stop the progress bar
        close(pb)
        parallel::stopCluster(cl)
      }
    } #end of GLMM method
    
    # Fit GEE (geeglm)
    if (method == 'gee') {
      if (quiet == FALSE) {
        sim.start <- Sys.time()
        geepack::geeglm(
          sim.dat[, 1] ~ trt,
          family = stats::binomial(link = 'logit'),
          id = clust,
          corstr = "exchangeable"
        )
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
        # initialize progress bar
        if (is.na(cores)) {
          prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta",
                                                 total = 5,
                                                 clear = FALSE)
          prog.bar$tick(0)
        }
      }
      require(foreach)
      if (!is.na(cores) & quiet == FALSE) {
        message("Fitting models")
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(2 / 5)
        Sys.sleep(1 / 100)
      }
      
      my.mod <- foreach::foreach(
        i = 1:nsim,
        .options.snow = opts,
        .packages = "geepack",
        .inorder = FALSE
      ) %fun% {
        geepack::geeglm(
          sim.dat[, i] ~ trt,
          family = stats::binomial(link = 'logit'),
          id = clust,
          corstr = "exchangeable"
        )
      }
      if (!is.na(cores) & quiet == FALSE) {
        message("Performing null model comparisons")
      }
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(3 / 5)
        Sys.sleep(1 / 100)
      }
      # get the overall p-values (>Chisq)
      model.compare <-
        foreach::foreach(i = 1:nsim, .inorder = FALSE) %fun% {
          anova(my.mod[[i]])
        }
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(4 / 5)
        Sys.sleep(1 / 100)
      }
      if (!is.na(cores) & quiet == FALSE) {
        message("Retrieving model summaries")
      }
      # get the model summaries
      model.values <-
        foreach::foreach(i = 1:nsim,
                         .packages = "car",
                         .inorder = FALSE) %fun% {
                           summary(my.mod[[i]])
                         }
      # turn off parallel computing
      if (!is.na(cores)) {
        #stop the progress bar
        close(pb)
        parallel::stopCluster(cl)
      }
    } # end of GEE method
    
    ## Output objects
    if (all.sim.data == TRUE) {
      complete.output.internal <-  list(
        "estimates" = model.values,
        "model.comparisons" = model.compare,
        "converged" = unlist(converged),
        "sim.data" = data.frame(trt, clust, sim.dat)
      )
    } else {
      complete.output.internal <-  list(
        "estimates" = model.values,
        "model.comparisons" = model.compare,
        "converged" = unlist(converged)
      )
    }
    return(complete.output.internal)
  } #end of function