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

test_that("plot function with extra conf int", {
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
  conf_int <- lm(formula = Postwt ~ Prewt + Treat + offset(Prewt),
                 data = anorexia) %>%
    tidy() %>%
    dplyr::mutate(lower = estimate - 1.96 * std.error,
                  upper = estimate + 1.96 * std.error) %>%
    dplyr::select(term, lower, upper) %>%
    data.frame()

  expect_no_error(plot(cs, extra_conf_int = conf_int))
})
