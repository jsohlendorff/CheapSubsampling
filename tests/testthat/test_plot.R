test_that("plot function", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  library(broom)
  fun <- function(data) {
    lm(formula = Postwt ~ Prewt + Treat + offset(Prewt), data = data) %>%
      tidy()
  }
  cs <- cheap_bootstrap(fun = fun, 
                        b = 20, 
                        data = anorexia,
                        est_col_name = "estimate",
                        par_col_names = "term")
  expect_no_error(plot(cs))
})