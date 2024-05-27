get_cheap_subsampling_ci <- function(est, boot_est, m_val, n_val, alpha) {
  b_val <- length(boot_est)
  s_val <- sqrt(mean((est - boot_est)^2))
  tq <- stats::qt(1 - alpha / 2, df = b_val)
  list(
    estimate = est,
    lower_b = est - tq * sqrt((m_val) / (n_val - m_val)) * s_val,
    upper_b = est + tq * sqrt((m_val) / (n_val - m_val)) * s_val
  )
}

##' Method implementing the cheap subsampling method for confidence intervals
##' 
##' Given a model object or a function that returns a vector of coefficients 
##' and a data set, 
##' this function computes confidence intervals using 
##' the cheap subsampling method.
##'
##' @title Cheap subsampling
##' @param x A function that returns a vector of coefficients 
##' or a model object which saves the call. In the second case, a
##' coef needs to be defined for 'class(x)'.
##' @param b Number of bootstrap samples.
##' @param m_val Subsample size. Defaults to 0.632 * nrow(data).
##' @param alpha Significance level. Defaults to 0.05.
##' @param parallelize Logical. 
##' Should the bootstrap samples be computed in parallel? 
##' Defaults to FALSE.
##' @param cores Number of cores to use for parallel computation. 
##' Defaults to detectCores().
##' @param data Data set to be used for the computation, if applicable. 
##' @return An object of class "cheap_subsampling" containing 
##' the point estimates and confidence intervals.
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs <- cheap_subsampling(x, b = 1000, data = anorexia)
##' cs
##' 
##' ## example with a model object
##' set.seed(123)
##' x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
##' cs2 <- cheap_subsampling(x, b = 1000)
##' cs2
##' 
##' ## example with a function call and parallel computation
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs3 <- cheap_subsampling(x, b = 1000, data = anorexia, parallelize = TRUE)
##' cs3
##' 
##' ## example with custom function returning coefficients 
##' ## example from ATE-package
##' library(survival)
##' library(riskRegression)
##' set.seed(10)
##' #### Survival settings  ####
##' #### ATE with Cox model ####
##' 
##' ## generate data
##' 
##' n <- 100
##' dtS <- sampleData(n, outcome="survival")
##' dtS$time <- round(dtS$time,1)
##' dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
##' 
##' ## estimate the Cox model
##' fit <- coxph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
##' ## fitter function which returns a named vector 
##' ate_fit_fun <- function(d) {
##'   ## main fitter function; ignore output
##'   invisible(capture.output(ate_fit <- summary(ate(fit, 
##'                                                   data = d, 
##'                                                   treatment = "X1", 
##'                                                   times = 5:8, 
##'                                                   se = FALSE))))
##'   ## extract the point estimates for risk difference
##'   res <- ate_fit$diffRisk$estimate
##'   ## name the point estimates
##'   names(res) <- paste0("ATE ",
##'                        ate_fit$diffRisk$A,
##'                        "-",
##'                        ate_fit$diffRisk$B, 
##'                        " (t = ",
##'                        ate_fit$diffRisk$time,")") 
##'   res
##' }
##' set.seed(102)
##' cs4 <- cheap_subsampling(ate_fit_fun, b = 100, data = dtS)
##' 
##' ## example from riskRegression with no existing coef function
##' set.seed(18)
##' learndat <- sampleData(200,outcome="binary")
##' testdat <- sampleData(1200,outcome="binary")
##'
##' ## score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
##' z<-Score(list("LR(X1+X2+X7+X9)"=lr1,
##'          "LR(X3+X5)"=lr2),
##'          formula=Y~1,
##'          data=testdat,
##'          metrics = "AUC")
##' coef.Score <- function(x) {
##'   res <- x$AUC$score$AUC
##'   names(res) <- paste0("AUC (model = ", x$AUC$score$model, ")")
##'   res
##' }
##' set.seed(102)
##' cs5 <- cheap_subsampling(z, b = 100)
##' }
cheap_subsampling <- function(x, b, m_val = NULL, alpha = 0.05, parallelize = FALSE, cores = parallel::detectCores(), data = NULL) {
  coef <- NULL
  if (!inherits(x, "function")) {
    ## tryCatch to retrieve x$call
    message("x is not a function, trying to retrieve call object from x")
    call <- x$call
    if (is.null(call)) {
      stop("x does not appear to have a call")
    }
    # get data from call in parent environment
    if (is.null(data)) {
      tryCatch({
        data <- eval(call$data, parent.frame())
      }, error = function(e) {
        stop("data not found in environment or missing from call object")
      })
    }
    ## make call object into string
    x <- function(d) {
      call$data <- d
      coef(eval(call))
    }
  }
  
  if (!inherits(data, "data.frame")) {
    stop("data needs to be a data frame")
  }
  n_val <- nrow(data)
  if (is.null(m_val)) {
    m_val <- round(0.632 * n_val)
  }
  est <- tryCatch({
    x(data)
  }, error = function(e) {
    stop("function x/derived from x failed with error: ", conditionMessage(e), ".\n Did you specify x as a function or as a model object with a corresponding call and coef function?")
  })

  ## check that est is a named vector
  if (!is.vector(est) || !is.numeric(est) || is.null(names(est))) {
    stop("function x/derived from x did not return a named vector of coefficients")
  }
  
  
  if (parallelize) {
    requireNamespace("parallel")
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("x", "b", "m_val", "alpha", "data"), envir = environment())
    tryCatch({
      boot_est <- parallel::parSapply(cl, seq_len(b), function(i) {
        x(data[sample(1:n_val, m_val, replace = FALSE), ])
      })
    }, error = function(e) {
      stop("bootstrap (parallel) computation failed with error: ", conditionMessage(e))
    })
    parallel::stopCluster(cl)
  } else {
    pb <- utils::txtProgressBar(min = 1, max = b, style = 3, width = 50, char = "=")
    tryCatch({
      boot_est <- sapply(seq_len(b), function(i) {
        utils::setTxtProgressBar(pb, i)
        x(data[sample(1:n_val, m_val, replace = FALSE), ])
      })
    }, error = function(e) {
      stop("bootstrap computation failed with error: ", conditionMessage(e))
    })
  }
  cat("\n")
  ## apply get_cheap_subsampling_ci for each row of boot_est, est
  tryCatch({
    res <- sapply(seq_len(length(est)), function(i) {
      get_cheap_subsampling_ci(est[i], boot_est[i, ], m_val, n_val, alpha)
    })
  }, error = function(e) {
    stop("computation of confidence intervals failed with error: ", conditionMessage(e), ". Does your function return a vector of coefficients?")
  })
  colnames(res) <- names(est)
  res <- list(res = res, boot_estimates = boot_est, b = b, m = m_val)
  class(res) <- "cheap_subsampling"
  res
}

##' Print method for cheap_subsampling objects
##' 
##' @title Print method for cheap_subsampling objects
##' @param x An object of class "cheap_subsampling"
##' @param ... Not applicable.
##' Prints the point estimates and confidence intervals.
##' @export
print.cheap_subsampling <- function(x, ...) {
  cat(paste0("Cheap subsampling results for subsample size m = ", x$m, " and ", x$b, " bootstrap samples\n"))
  print(x$res, ...)
  invisible(x)
}
