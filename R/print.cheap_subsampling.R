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
##' Print method for cheap_bootstrap objects
##' 
##' @title Print method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Prints the point estimates and confidence intervals.
##' @export
print.cheap_bootstrap <- function(x, ...) {
  if (x$type == "subsampling") {
    cat(paste0("Cheap subsampling results for subsample size m = ", x$size, " and ", x$b, " bootstrap samples\n"))
  } else {
    cat(paste0("Cheap (non-parametric) bootstrap results for ", x$b, " bootstrap samples\n"))
  }
  print(x$res, ...)
  invisible(x)
}

######################################################################
### print.cheap_subsampling.R ends here
