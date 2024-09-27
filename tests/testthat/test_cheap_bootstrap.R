## make unit tests
test_that("cheap subsampling with function works", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  }
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, data = anorexia)
  expect_output(print(cs))
  expect_output(summary(cs))
  expect_output(print(cs$res))
})

test_that("cheap subsampling with model work", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10)
  expect_output(print(cs))
})

test_that("cheap subsampling in parallel runs (function)", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  }
  set.seed(6)
  cs <- cheap_bootstrap(
    x,
    b = 10,
    data = anorexia,
    parallel_args = list(
      parallelize = TRUE,
      cores = 2
    )
  )
  expect_output(print(cs))
})

test_that("cheap subsampling in parallel runs (model)", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x,
    b = 10,
    parallel_args = list(parallelize = TRUE, cores = 2)
  )
  expect_output(print(cs))
})

test_that("cheap bootstrap with function works", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  }
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, data = anorexia, type = "non_parametric")
  expect_output(print(cs))
})

test_that("cheap bootstrap with model works", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, type = "non_parametric")
  expect_output(print(cs))
})

test_that("error when x is not a function, but doesn't contain a call", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  x$call <- NULL
  set.seed(6)
  expect_error(
    cheap_bootstrap(x, b = 10),
    "fun does not appear to have a call"
  )
})

test_that("error when coef function is not defined", {
  require(riskRegression)
  set.seed(18)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")

  ## score logistic regression models
  lr1 <- glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 <- glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  set.seed(102)
  expect_error(
    cheap_bootstrap(z, b = 10),
    "function fun/derived from fun did not return a named vector of coefficients"
  )
})

test_that("error when data does not exist in parent environment", {
  set.seed(18)
  require(riskRegression)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")

  ## score logistic regression models
  lr1 <- glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 <- glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  remove(testdat)
  set.seed(102)
  expect_error(
    cheap_bootstrap(z, b = 10),
    "data not found in environment or missing from call object"
  )
})

test_that("error when data is not in correct format", {
  set.seed(18)
  require(riskRegression)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")

  ## score logistic regression models
  lr1 <- glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 <- glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  testdat <- 2
  set.seed(102)
  expect_error(
    cheap_bootstrap(z, b = 10),
    "Data needs to be a data frame"
  )
  expect_error(
    get_data_and_call(z, data = NULL, environment()),
    "Data needs to be a data frame"
  )
})

test_that("error when coef gives error", {
  set.seed(18)
  require(riskRegression)
  learndat <- sampleData(200, outcome = "binary")
  testdat <- sampleData(1200, outcome = "binary")

  ## score logistic regression models
  lr1 <- glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 <- glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  coef.Score <- function(x) {
    stop("oh no")
  }
  set.seed(102)
  expect_error(
    cheap_bootstrap(z, b = 10),
    "function fun/derived from fun failed with error: oh no.\n Did you specify fun as a function or as a model object with a corresponding call and coef function?"
  )
})

test_that("error when coef does not return a named vector", {
  set.seed(18)
  require(riskRegression)
  learndat <- sampleData(200, outcome = "binary")
  testdat <- sampleData(1200, outcome = "binary")

  ## score logistic regression models
  lr1 <- glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 <- glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  coef.Score <- function(x) {
    list(a = 2)
  }
  set.seed(102)
  expect_error(
    cheap_bootstrap(z, b = 10),
    "function fun/derived from fun did not return a named vector of coefficients"
  )
})

test_that("error when bootstrap gives error", {
  require(survival)
  require(riskRegression)
  #### Survival settings  ####
  #### ATE with Cox model ####

  ## generate data
  set.seed(16)
  n <- 100
  dt <- sampleData(n, outcome = "survival")
  dt$time <- round(dt$time, 1)
  dt$X1 <- factor(
    rbinom(n,
      prob = c(0.3, 0.4),
      size = 2
    ),
    labels = paste0("T", 0:2)
  )

  ## estimate the Cox model
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    fit <- coxph(
      formula = Surv(time, event) ~ X1 + X2,
      data = d,
      y = TRUE,
      x = TRUE
    )
    ## main fitter function; ignore output
    invisible(capture.output(suppressWarnings(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    ))))
    ## extract the point estimates for risk difference
    res <- ate_fit$diffRisk$estimate
    ## name the point estimates
    names(res) <- paste0(
      "ATE ",
      ate_fit$diffRisk$A,
      "-",
      ate_fit$diffRisk$B,
      " (t = ",
      ate_fit$diffRisk$time,
      ")"
    )
    if (runif(1) > 0.5) {
      stop("oh no!")
    }
    res
  }
  set.seed(105)
  expect_error(
    cheap_bootstrap(ate_fit_fun, b = 10, data = dt),
    "Bootstrap computation failed with error: oh no!"
  )
})

test_that("error when bootstrap gives error (parallel)", {
  require(survival)
  require(riskRegression)
  #### Survival settings  ####
  #### ATE with Cox model ####

  ## generate data
  set.seed(16)
  n <- 100
  dt <- sampleData(n, outcome = "survival")
  dt$time <- round(dt$time, 1)
  dt$X1 <- factor(
    rbinom(n,
      prob = c(0.3, 0.4),
      size = 2
    ),
    labels = paste0("T", 0:2)
  )

  ## estimate the Cox model
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    require(survival)
    require(riskRegression)
    fit <- coxph(
      formula = Surv(time, event) ~ X1 + X2,
      data = d,
      y = TRUE,
      x = TRUE
    )
    ## main fitter function; ignore output
    invisible(capture.output(suppressWarnings(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    ))))
    ## extract the point estimates for risk difference
    res <- ate_fit$diffRisk$estimate
    ## name the point estimates
    names(res) <- paste0(
      "ATE ",
      ate_fit$diffRisk$A,
      "-",
      ate_fit$diffRisk$B,
      " (t = ",
      ate_fit$diffRisk$time,
      ")"
    )
    if (runif(1) > 0.5) {
      stop("oh no!")
    }
    res
  }
  set.seed(105)
  expect_error(
    suppressWarnings(cheap_bootstrap(ate_fit_fun,
      b = 10,
      data = dt,
      parallel_args = list(
        parallelize = TRUE,
        cores = 2
      )
    ))
  )
})

test_that("more than one argument for parallelize", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  }
  parallel_args <- list(parallelize = c(FALSE, TRUE))
  expect_error(
    cheap_bootstrap(x,
      parallel_args = parallel_args,
      data = anorexia
    ),
    "parallelize must be of length 1"
  )
})

test_that("calling fun fails", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    stop("oh no")
  }
  expect_error(
    cheap_bootstrap(x, data = anorexia),
    "function fun/derived from fun failed with error: oh no.\n Did you specify fun as a function or as a model object with a corresponding call and coef function?"
  )
})

test_that("fail in get_cheap_subsampling_confidence_interval", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    if (runif(1) < 0.5) {
      c(a = 1, b = 2, b = 3)
    } else {
      list(a = 1, b = 2)
    }
  }
  set.seed(8)
  expect_error(
    suppressWarnings(cheap_bootstrap(x, data = anorexia)),
    "Computation of confidence intervals failed with error: non-numeric argument to binary operator. Does your function return a vector of coefficients?"
  )
})

test_that("adaptive cheap_bootstrap functionality", {
  ## example with adaptive method
  utils::data(anorexia, package = "MASS")
  set.seed(123)
  fun <- function(data) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = data))
  }
  cs <- cheap_bootstrap(
    fun = fun,
    data = anorexia,
    adapt_args = list(
      adapt = TRUE,
      precision = 0.1,
      max_b = 200
    )
  )
  expect_output(print(cs))
})

test_that("error when bootstrap gives error (adapt)", {
  require(survival)
  require(riskRegression)
  #### Survival settings  ####
  #### ATE with Cox model ####

  ## generate data
  set.seed(16)
  n <- 100
  dt <- sampleData(n, outcome = "survival")
  dt$time <- round(dt$time, 1)
  dt$X1 <- factor(
    rbinom(n,
      prob = c(0.3, 0.4),
      size = 2
    ),
    labels = paste0("T", 0:2)
  )

  ## estimate the Cox model
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    fit <- coxph(
      formula = Surv(time, event) ~ X1 + X2,
      data = d,
      y = TRUE,
      x = TRUE
    )
    ## main fitter function; ignore output
    invisible(capture.output(suppressWarnings(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    ))))
    ## extract the point estimates for risk difference
    res <- ate_fit$diffRisk$estimate
    ## name the point estimates
    names(res) <- paste0(
      "ATE ",
      ate_fit$diffRisk$A,
      "-",
      ate_fit$diffRisk$B,
      " (t = ",
      ate_fit$diffRisk$time,
      ")"
    )
    if (runif(1) > 0.5) {
      stop("oh no!")
    }
    res
  }
  set.seed(105)
  expect_error(
    cheap_bootstrap(ate_fit_fun,
      adapt_args = list(
        adapt = TRUE,
        precision = 0.1,
        max_b = 200
      ),
      data = dt
    ),
    "Bootstrap computation failed with error: oh no! for bootstrap iteration b = 2"
  )
})


test_that("fail in get_cheap_subsampling_confidence_interval (adapt)", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    if (runif(1) < 0.5) {
      c(a = 1, b = 2, b = 3)
    } else {
      c(a = 1, b = 2)
    }
  }
  set.seed(12)
  expect_error(suppressWarnings(cheap_bootstrap(x, data = anorexia, adapt = TRUE)))
})

test_that("get all bootstrap interval", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  cs <- cheap_bootstrap(x, b = 10, data = anorexia)
  cs <- get_bootstrap_cis(cs)
  expect_output(print(cs))
})

test_that("get all bootstrap interval", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  cs <- cheap_bootstrap(x, b = 10, data = anorexia)
  cs$boot_estimates <- list(a = 1, b = 2)
  expect_error(
    get_bootstrap_cis(cs),
    "Computation of confidence intervals failed with error: incorrect number of dimensions for b = 1"
  )
})

test_that("get all bootstrap interval", {
  utils::data(anorexia, package = "MASS")
  ## example with a function call
  set.seed(123)
  x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  cs <- cheap_bootstrap(x,
    b = 10,
    data = anorexia,
    additional_args = list(keep_boot_estimates = FALSE, verbose = FALSE)
  )
  expect_error(get_bootstrap_cis(cs), "No bootstrap estimates found in x")
})

test_that("fail in get_cheap_subsampling_confidence_interval (adapt)", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    if (runif(1) < 0.5) {
      c(a = 1, b = 2, b = 3)
    } else {
      list(a = 1, b = 2)
    }
  }
  set.seed(8)
  ## note does not let me test
  expect_error(
    suppressWarnings(cheap_bootstrap(x,
      data = anorexia,
      adapt_args = list(adapt = TRUE, precision = 0.1, max_b = 200)
    ))
  )
})
