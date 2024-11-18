test_that("plot function, with extra_conf_int, null extra_conf_int", {
  utils::data(anorexia, package = "MASS")
  set.seed(123)
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
  conf_int <- lm(formula = Postwt ~ Prewt + Treat + offset(Prewt),
                 data = anorexia) %>%
    tidy() %>%
    dplyr::mutate(lower = estimate - 1.96 * std.error,
                  upper = estimate + 1.96 * std.error) %>%
    dplyr::select(term, lower, upper) %>%
    data.frame()
  ## example additional confidence interval
  expect_no_error(plot(cs, extra_conf_int = conf_int))
  
  ## examples with additional confidence interval incorrectly specified
  expect_error(plot(cs, extra_conf_int = list()))
  expect_error(plot(cs, extra_conf_int = data.frame(x=3)))
})
