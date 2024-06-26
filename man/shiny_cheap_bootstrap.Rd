% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{shiny_cheap_bootstrap}
\alias{shiny_cheap_bootstrap}
\title{Shiny app for cheap subsampling}
\usage{
shiny_cheap_bootstrap(fun, data = NULL)
}
\arguments{
\item{fun}{A function that returns a vector of coefficients
or a model object which saves the call. In the second case, a
coef needs to be defined for 'class(fun)'.}

\item{data}{Data set to be used for the computation, if applicable.}
}
\description{
Shiny app for cheap subsampling
}
\details{
This function creates a shiny app for the cheap subsampling method.
The user can vary the number of bootstrap samples and the subsample size.
}
\examples{
\dontrun{
utils::data(anorexia, package = "MASS")
## example with a function call
set.seed(123)
fun <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt),
                           data = d))
shiny_cheap_bootstrap(fun, data = anorexia)
}
}
