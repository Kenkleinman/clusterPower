#' Look up which internal functions are called by exported functions.
#'
#' @author Alexandria Sakrejda (\email{acbro0@@umass.edu})
#'
#' @param packageName The name of the package in quotes. Defaults to "clusterPower".
#'
#' @return List of internal functions and the line numbers in which they appear 
#' inside external functions.
#'
#' @export

package_map_helper <-
  function(packageName = "clusterPower") {
    toSearchWithin <-
      function(exported_fxn_name) {
        temp <- deparse(utils::getAnywhere(exported_fxn_name)[[2]][[1]])
        names(temp) <- exported_fxn_name
        return(temp)
      }
    check <- function(x, pattern) {
      out <- NULL
      if (any(grepl(x = x, pattern = pattern) == TRUE)) {
        out <- grep(x = x, pattern = pattern)
      }
      return(out)
    }
    
    organize <-
      function(pattern = pattern,
               toSearch = holder,
               names = y) {
        argM <- lapply(toSearch, check, pattern = pattern)
        names(argM) <- names
        argM <- unlist(argM)
        return(argM)
      }
    
    findFortran <- function(allx) {
      out <-
        any(grepl(x = attributes(utils::getAnywhere(allx)[["objs"]][[1]])$class, pattern = "FortranRoutine"))
      return(out)
    }
    
    pck <- paste0("package:", packageName, sep = "")
    y <- utils::lsf.str(pck)
    attributes(y) <- NULL
    holder <- list()
    holder <- lapply(y, toSearchWithin)
    all <- ls(getNamespace(packageName), all.names = FALSE)
    fortranvec <- unlist(lapply(all, findFortran))
    all <- all[fortranvec == FALSE]
    toSearchFor <- setdiff(all, y)
    q <- lapply(toSearchFor, organize)
    names(q) <- toSearchFor
    return(q)
  }
