test_that("test cheap subsampling", {
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
                        par_col_names = "term",
                        size = round(nrow(anorexia) * 0.8))
  expect_no_error(summary(cs))
})

test_that("test cheap bootstrap", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  library(broom)
  fun <- function(data) {
    lm(formula = Postwt ~ Prewt + Treat + offset(Prewt), data = data) %>%
      tidy()
  }
  cb <- cheap_bootstrap(fun = fun, 
                        b = 20, 
                        data = anorexia,
                        est_col_name = "estimate",
                        par_col_names = "term",
                        type = "non_parametric")
  expect_no_error(summary(cb))
})