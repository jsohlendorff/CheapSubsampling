test_that("plot function", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  cs <- cheap_bootstrap(x, b = 1000, data = anorexia)
  expect_no_error(plot(cs))
})

test_that("plot function (bootstrap)", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  cs <- cheap_bootstrap(x, b = 1000, data = anorexia, type = "non_parametric")
  expect_no_error(plot(cs))
})