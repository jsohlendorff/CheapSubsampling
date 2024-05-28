#' CheapSubsampling: Fast Subsampling
#'
#' @keywords internal
# "_PACKAGE"
#' @name CheapSubsampling-package
#' @aliases CheapSubsampling-package
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel detectCores makeCluster parSapply stopCluster
#' @importFrom shiny checkboxInput fluidPage mainPanel renderTable selectInput sidebarLayout sliderInput tableOutput titlePanel observeEvent
#' @importFrom ggplot2 ggplot geom_point geom_line geom_smooth theme_minimal aes
NULL

# useful commands
# devtools::create("CheapSubsampling")
# usethis::use_gpl_license(version = 2) ## to add license
# devtools::test() ## to make tests and test
# usethis::use_github_action() ## to add github actions
# covr::package_coverage() ## to check coverage
# zero_coverage(package_coverage()) ## to check which lines are not covered
# usethis::use_coverage()
# lintr?
# note: make shiny app