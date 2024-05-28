## make unit tests
test_that("cheap_subsampling_function", {
  utils::data(anorexia, package = "MASS")
  x <- function(d)
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, data = anorexia)
  check <- structure(
    list(
      c(`(Intercept)` = 49.7711090149846),
      c(`(Intercept)` = 15.4837723109149),
      c(`(Intercept)` = 84.0584457190543),
      c(Prewt = -0.565538849639096),
      c(Prewt = -0.964312324218352),
      c(Prewt = -0.16676537505984),
      c(TreatCont = -4.0970655280729),
      c(TreatCont = -8.50522658403253),
      c(TreatCont = 0.31109552788673),
      c(TreatFT = 4.56306265291879),
      c(TreatFT = -0.420345452050825),
      c(TreatFT = 9.5464707578884)
    ),
    dim = 3:4,
    dimnames = list(
      c("estimate", "cheap_lower", "cheap_upper"),
      c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
    )
  )
  expect_equal(cs$res, check)
})

test_that("cheap_subsampling_model", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10)
  check <- structure(
    list(
      c(`(Intercept)` = 49.7711090149846),
      c(`(Intercept)` = 15.4837723109149),
      c(`(Intercept)` = 84.0584457190543),
      c(Prewt = -0.565538849639096),
      c(Prewt = -0.964312324218352),
      c(Prewt = -0.16676537505984),
      c(TreatCont = -4.0970655280729),
      c(TreatCont = -8.50522658403253),
      c(TreatCont = 0.31109552788673),
      c(TreatFT = 4.56306265291879),
      c(TreatFT = -0.420345452050825),
      c(TreatFT = 9.5464707578884)
    ),
    dim = 3:4,
    dimnames = list(
      c("estimate", "cheap_lower", "cheap_upper"),
      c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
    )
  )
  expect_equal(cs$res, check)
})

test_that("cheap_subsampling_function_parallel", {
  utils::data(anorexia, package = "MASS")
  x <- function(d)
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  set.seed(6)
  cs <- cheap_bootstrap(
    x,
    b = 10,
    data = anorexia,
    parallelize = TRUE,
    cores = 2
  )
  expect_output(print(cs))
})

test_that("cheap_subsampling_model_parallel", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x,
                        b = 10,
                        parallelize = TRUE,
                        cores = 2)
  expect_output(print(cs))
})

test_that("cheap_bootstrap_function", {
  utils::data(anorexia, package = "MASS")
  x <- function(d)
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, data = anorexia, type = "non_parametric")
  check <- structure(
    list(
      c(`(Intercept)` = 49.7711090149846),
      c(`(Intercept)` = 9.04028074218562),
      c(`(Intercept)` = 90.5019372877836),
      c(Prewt = -0.565538849639096),
      c(Prewt = -1.06914141201443),
      c(Prewt = -0.0619362872637669),
      c(TreatCont = -4.0970655280729),
      c(TreatCont = -7.86760563700546),
      c(TreatCont = -0.326525419140337),
      c(TreatFT = 4.56306265291879),
      c(TreatFT = -1.12128954094085),
      c(TreatFT = 10.2474148467784)
    ),
    dim = 3:4,
    dimnames = list(
      c("estimate", "cheap_lower", "cheap_upper"),
      c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
    )
  )
  expect_equal(cs$res, check)
  expect_output(print(cs))
})


test_that("error when x is not a function, but doesn't contain a call", {
  require(survival)
  require(riskRegression)
  set.seed(10)
  #### Survival settings  ####
  #### ATE with Cox model ####
  
  ## generate data
  
  n <- 100
  dtS <- sampleData(n, outcome = "survival")
  dtS$time <- round(dtS$time, 1)
  dtS$X1 <- factor(rbinom(n, prob = c(0.3, 0.4) , size = 2), labels = paste0("T", 0:2))
  
  ## estimate the Cox model
  fit <- coxph(
    formula = Surv(time, event) ~ X1 + X2,
    data = dtS,
    y = TRUE,
    x = TRUE
  )
  ate_fit <- ate(
    fit,
    data = dtS,
    treatment = "X1",
    times = 5:8,
    se = FALSE,
    verbose = FALSE
  )
  expect_error(cheap_bootstrap(ate_fit, b = 10))
})

test_that("error when coef function is not defined", {
  set.seed(18)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")
  
  ## score logistic regression models
  lr1 = glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 = glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  set.seed(102)
  expect_error(cheap_bootstrap(z, b = 10))
})

test_that("error when data does not exist in parent environment", {
  set.seed(18)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")
  
  ## score logistic regression models
  lr1 = glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 = glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  remove(testdat)
  set.seed(102)
  expect_error(cheap_bootstrap(z, b = 10))
})

test_that("error when data is not in correct format", {
  set.seed(18)
  learndat <- sampleData(20, outcome = "binary")
  testdat <- sampleData(120, outcome = "binary")
  
  ## score logistic regression models
  lr1 = glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 = glm(Y ~ X3 + X5, data = learndat, family = binomial)
  z <- Score(
    list("LR(X1+X2+X7+X9)" = lr1, "LR(X3+X5)" = lr2),
    formula = Y ~ 1,
    data = testdat,
    metrics = "AUC"
  )
  testdat <- 2
  set.seed(102)
  expect_error(cheap_bootstrap(z, b = 10))
})

test_that("error when coef gives error", {
  set.seed(18)
  learndat <- sampleData(200, outcome = "binary")
  testdat <- sampleData(1200, outcome = "binary")
  
  ## score logistic regression models
  lr1 = glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 = glm(Y ~ X3 + X5, data = learndat, family = binomial)
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
  expect_error(cheap_bootstrap(z, b = 10))
})

test_that("error when coef does not return a named vector", {
  set.seed(18)
  learndat <- sampleData(200, outcome = "binary")
  testdat <- sampleData(1200, outcome = "binary")
  
  ## score logistic regression models
  lr1 = glm(Y ~ X1 + X2 + X7 + X9, data = learndat, family = binomial)
  lr2 = glm(Y ~ X3 + X5, data = learndat, family = binomial)
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
  expect_error(cheap_bootstrap(z, b = 10))
})

test_that("error when bootstrap gives weird results", {
  require(survival)
  require(riskRegression)
  set.seed(10)
  #### Survival settings  ####
  #### ATE with Cox model ####
  
  ## generate data
  
  n <- 100
  dtS <- sampleData(n, outcome = "survival")
  dtS$time <- round(dtS$time, 1)
  dtS$X1 <- factor(rbinom(n, prob = c(0.3, 0.4) , size = 2), labels = paste0("T", 0:2))
  
  ## estimate the Cox model
  fit <- coxph(
    formula = Surv(time, event) ~ X1 + X2,
    data = dtS,
    y = TRUE,
    x = TRUE
  )
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    ## main fitter function; ignore output
    invisible(capture.output(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    )))
    res <- ate_fit$diffRisk$estimate ## extract the point estimates for risk difference
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
      res <- list(a = 2)
    }
    res
  }
  set.seed(102)
  expect_error(cheap_bootstrap(ate_fit_fun, b = 10, data = dtS))
})

test_that("error when bootstrap gives error", {
  require(survival)
  require(riskRegression)
  set.seed(10)
  #### Survival settings  ####
  #### ATE with Cox model ####
  
  ## generate data
  
  n <- 100
  dtS <- sampleData(n, outcome = "survival")
  dtS$time <- round(dtS$time, 1)
  dtS$X1 <- factor(rbinom(n, prob = c(0.3, 0.4) , size = 2), labels = paste0("T", 0:2))
  
  ## estimate the Cox model
  fit <- coxph(
    formula = Surv(time, event) ~ X1 + X2,
    data = dtS,
    y = TRUE,
    x = TRUE
  )
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    ## main fitter function; ignore output
    invisible(capture.output(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    )))
    res <- ate_fit$diffRisk$estimate ## extract the point estimates for risk difference
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
  set.seed(102)
  expect_error(cheap_bootstrap(ate_fit_fun, b = 10, data = dtS))
})

test_that("error when bootstrap gives error (parallel)", {
  require(survival)
  require(riskRegression)
  set.seed(10)
  #### Survival settings  ####
  #### ATE with Cox model ####
  
  ## generate data
  
  n <- 100
  dtS <- sampleData(n, outcome = "survival")
  dtS$time <- round(dtS$time, 1)
  dtS$X1 <- factor(rbinom(n, prob = c(0.3, 0.4) , size = 2), labels = paste0("T", 0:2))
  
  ## estimate the Cox model
  fit <- coxph(
    formula = Surv(time, event) ~ X1 + X2,
    data = dtS,
    y = TRUE,
    x = TRUE
  )
  ## fitter function which returns a named vector
  ate_fit_fun <- function(d) {
    ## main fitter function; ignore output
    invisible(capture.output(ate_fit <- summary(
      ate(
        fit,
        data = d,
        treatment = "X1",
        times = 5:8,
        se = FALSE
      )
    )))
    res <- ate_fit$diffRisk$estimate ## extract the point estimates for risk difference
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
  set.seed(102)
  expect_error(cheap_bootstrap(
    ate_fit_fun,
    b = 10,
    data = dtS,
    parallelize = TRUE,
    cores = 2
  ))
})