test_that("table function shiny", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  expect_no_error(shiny_cheap_bootstrap(x, data = anorexia))
})

test_that("plot function shiny", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  expect_no_error(shiny_cheap_bootstrap_plot(x, data = anorexia))
})