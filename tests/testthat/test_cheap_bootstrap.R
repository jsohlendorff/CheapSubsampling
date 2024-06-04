## make unit tests
test_that("cheap subsampling with function works", {
  utils::data(anorexia, package = "MASS")
  x <- function(d) {
    coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
  }
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, data = anorexia)
  check <- data.frame(
    estimate = c(
      49.7711090149846, -0.565538849639096,
      -4.0970655280729,
      4.56306265291879
    ),
    cheap_lower = c(
      15.4837723109149,
      -0.964312324218352, -8.50522658403253, -0.420345452050825
    ),
    cheap_upper = c(
      84.0584457190543,
      -0.16676537505984,
      0.31109552788673,
      9.5464707578884
    ),
    parameter = c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
  )
  expect_output(print(cs))
  expect_output(summary(cs))
  expect_equal(cs$res, check)
})

test_that("cheap subsampling with model work", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10)
  check <- data.frame(
    estimate = c(
      49.7711090149846, -0.565538849639096,
      -4.0970655280729,
      4.56306265291879
    ),
    cheap_lower = c(
      15.4837723109149,
      -0.964312324218352, -8.50522658403253, -0.420345452050825
    ),
    cheap_upper = c(
      84.0584457190543,
      -0.16676537505984,
      0.31109552788673,
      9.5464707578884
    ),
    parameter = c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
  )
  expect_equal(cs$res, check)
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
  check <- data.frame(
    estimate = c(
      49.7711090149846,
      -0.565538849639096, -4.0970655280729,
      4.56306265291879
    ),
    cheap_lower = c(
      9.04028074218562, -1.06914141201443,
      -7.86760563700546,
      -1.12128954094085
    ),
    cheap_upper = c(
      90.5019372877836, -0.0619362872637669,
      -0.326525419140337,
      10.2474148467784
    ),
    parameter = c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
  )
  expect_equal(cs$res, check)
  expect_output(print(cs))
})

test_that("cheap bootstrap with model works", {
  utils::data(anorexia, package = "MASS")
  x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
  set.seed(6)
  cs <- cheap_bootstrap(x, b = 10, type = "non_parametric")
  check <- data.frame(
    estimate = c(
      49.7711090149846,
      -0.565538849639096, -4.0970655280729,
      4.56306265291879
    ),
    cheap_lower = c(
      9.04028074218562, -1.06914141201443,
      -7.86760563700546,
      -1.12128954094085
    ),
    cheap_upper = c(
      90.5019372877836, -0.0619362872637669,
      -0.326525419140337,
      10.2474148467784
    ),
    parameter = c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
  )
  expect_equal(cs$res, check)
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
    cheap_bootstrap(ate_fit_fun,
      b = 10,
      data = dt,
      parallel_args = list(
        parallelize = TRUE,
        cores = 2
      )
    ),
    "Bootstrap computation failed with error: 2 nodes produced errors; first error: oh no!"
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
  check <- data.frame(
    estimate = c(
      49.7711090149846,
      -0.565538849639096, -4.0970655280729,
      4.56306265291879
    ),
    cheap_lower = c(
      14.6373754108099, -0.995226714362994,
      -7.63484151214865,
      -0.336391031355021
    ),
    cheap_upper = c(
      84.9048426191593, -0.135850984915198,
      -0.559289543997145,
      9.46251633719259
    ),
    parameter = c("(Intercept)", "Prewt", "TreatCont", "TreatFT")
  )
  expect_equal(cs$res, check)
  set.seed(5)
  cs <- cheap_bootstrap(
    fun = fun,
    data = anorexia,
    adapt_args = list(
      adapt = TRUE,
      precision = 0.001,
      max_b = 20
    )
  )
  set.seed(5)
  cs2 <- cheap_bootstrap(
    fun = fun,
    data = anorexia,
    b = 20,
  )
  expect_equal(cs$res, cs2$res)
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
  check <- structure(
    list(
      estimate = c(
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879,
        49.7711090149846,
        -0.565538849639096, -4.0970655280729,
        4.56306265291879
      ),
      cheap_lower = c(
        -79.4285481289019, -2.14168407474334,
        -6.83259459748369,
        -24.9105778836177,
        -17.1349397365194, -1.38492462558751,
        -9.58065183489002,
        -9.8094121820358,
        -7.99606974260923, -1.26214389858383,
        -7.86034919029748,
        -6.81305174306396,
        -0.373611667614128, -1.16460370513152,
        -7.17502349457684,
        -4.11635652813816,
        2.08264275995942, -1.15415769253428,
        -7.60199813788498,
        -4.42139768390172,
        7.81489286531418, -1.08165408018403,
        -7.15910098356241,
        -3.2907927362683,
        11.4447469284138, -1.03547780156731,
        -6.95982942370666,
        -2.47244596284023,
        12.5040077822593, -1.01608881646726,
        -6.83584860815544,
        -2.42113556676505,
        13.1030560803521, -1.00428979090368,
        -6.68033667893457,
        -2.15833977608044,
        4.29275730868622, -1.11977636811514,
        -6.56115659000526,
        -2.09675989574963
      ),
      cheap_upper = c(
        178.970766158871,
        1.01060637546515,
        -1.3615364586621,
        34.0367031894552,
        116.677157766489,
        0.253846926309322,
        1.38652077874422,
        18.9355374878734,
        107.538287772578,
        0.131066199305641,
        -0.333781865848315,
        15.9391770489015,
        99.9158296975833,
        0.0335260058533257,
        -1.01910756156896,
        13.2424818339757,
        97.4595752700098,
        0.0230799932560849,
        -0.592132918260816,
        13.5475229897393,
        91.727325164655, -0.0494236190941673,
        -1.03503007258339,
        12.4169180421059,
        88.0974711015554, -0.0955998977108781,
        -1.23430163243914,
        11.5985712686778,
        87.0382102477099, -0.114988882810934,
        -1.35828244799036,
        11.5472608726026,
        86.4391619496171, -0.126787908374515,
        -1.51379437721123,
        11.284465081918,
        95.249460721283, -0.011301331163057,
        -1.63297446614053,
        11.2228852015872
      ),
      parameter = c(
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT",
        "(Intercept)",
        "Prewt",
        "TreatCont",
        "TreatFT"
      ),
      b = c(
        1L,
        1L,
        1L,
        1L,
        2L,
        2L,
        2L,
        2L,
        3L,
        3L,
        3L,
        3L,
        4L,
        4L,
        4L,
        4L,
        5L,
        5L,
        5L,
        5L,
        6L,
        6L,
        6L,
        6L,
        7L,
        7L,
        7L,
        7L,
        8L,
        8L,
        8L,
        8L,
        9L,
        9L,
        9L,
        9L,
        10L,
        10L,
        10L,
        10L
      )
    ),
    row.names = c(NA, -40L),
    class = "data.frame"
  )
  expect_equal(cs, check)
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
