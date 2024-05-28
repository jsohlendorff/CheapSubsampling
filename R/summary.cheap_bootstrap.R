### summary.cheap_bootstrap.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 28 2024 (08:39)
## Version:
## Last-Updated: May 28 2024 (13:20)
##           By: Thomas Alexander Gerds
##     Update #: 4
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
##' @param object An object of class "cheap_subsampling"
##' @param ... Not applicable.
##' @return Summary table with the point estimates and confidence intervals.
##' @export

summary.cheap_bootstrap <- function(object,print = TRUE, ...) {
  y <- object$res
  if (print[[1]] == TRUE){
    print(y,...)
  }
  invisible(y)
}

######################################################################
### summary.cheap_subsampling.R ends here
