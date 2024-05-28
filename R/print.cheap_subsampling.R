### print.cheap_subsampling.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 28 2024 (08:39) 
## Version: 
## Last-Updated: May 28 2024 (08:39) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Print method for cheap_subsampling objects
##' 
##' @title Print method for cheap_subsampling objects
##' @param x An object of class "cheap_subsampling"
##' @param ... Not applicable.
##' Prints the point estimates and confidence intervals.
##' @export
print.cheap_subsampling <- function(x, ...) {
  cat(paste0("Cheap subsampling results for subsample size m = ", x$m, " and ", x$b, " bootstrap samples\n"))
  print(x$res, ...)
  invisible(x)
}


######################################################################
### print.cheap_subsampling.R ends here
