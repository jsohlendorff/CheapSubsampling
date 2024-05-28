### summary.cheap_subsampling.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 28 2024 (08:39) 
## Version: 
## Last-Updated: May 28 2024 (12:57) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Summary method for cheap_subsampling objects
##' 
##' @title Summary method for cheap_subsampling objects
##' @param x An object of class "cheap_subsampling"
##' @param ... Not applicable.
##' @return Summary table with the point estimates and confidence intervals.
##' @export
summary.cheap_subsampling <- function(x, ...) {
  cat(paste0("Cheap subsampling results for subsample size m = ", x$m, " and ", x$b, " bootstrap samples\n"))
  summary(x$res, ...)
  invisible(x)
}


######################################################################
### summary.cheap_subsampling.R ends here
