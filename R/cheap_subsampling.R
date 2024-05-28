##' Method implementing the cheap subsampling method for confidence intervals
##' 
##' Given a model object or a function that returns a vector of coefficients 
##' and a data set, 
##' this function computes confidence intervals using 
##' the cheap subsampling method.
##'
##' @title Cheap subsampling
##' @param fun
##' @param b Number of bootstrap samples.
##' @param size Subsample size. Defaults to 0.632 * nrow(data).
##' @param alpha Significance level. Defaults to 0.05.
##' @param keep_boot_estimates Logical. If \code{TRUE} add the subsampling bootstrap estimates to the output list. 
##' @param parallelize Logical.  Should the bootstrap samples be
##'     computed in parallel?  Defaults to FALSE.
##' @param cores Number of cores to use for parallel computation.
##'     Defaults to detectCores().
##' @param data Data set to be used for the computation, if
##'     applicable.
##' @param x A function that returns a vector of coefficients or a
##'     model object which saves the call. In the second case, a coef
##'     needs to be defined for 'class(x)'.
##' @return An object of class "cheap_subsampling" containing the
##'     point estimates and confidence intervals.
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' fun <- function(data) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = data))
##' cs <- cheap_subsampling(fun=fun, b = 1000, data = anorexia)
##' cs
##' 
##' ## example with a model object
##' set.seed(123)
##' x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
##' cs2 <- cheap_subsampling(fun=x, b = 1000)
##' cs2
##' 
##' ## example with a function call and parallel computation
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs3 <- cheap_subsampling(fun=x, b = 1000, data = anorexia, parallelize = TRUE)
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
##' cs4
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
cheap_subsampling <- function(fun,
                              b,
                              size = NULL,
                              alpha = 0.05,
                              keep_boot_estimates = FALSE,
                              parallelize = FALSE,
                              cores = parallel::detectCores(),
                              data = NULL) {
    coef <- NULL
    if (!inherits(fun, "function")) {
        ## tryCatch to retrieve fun$call
        message("fun is not a function, trying to retrieve call object from fun")
        call <- fun$call
        if (is.null(call)) {
            stop("fun does not appear to have a call")
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
        fun <- function(d) {
            call$data <- d
            coef(eval(call))
        }
    }
    
    if (!inherits(data, "data.frame")) {
        stop("data needs to be a data frame")
    }
    n_val <- nrow(data)
    if (is.null(size)) {
        size <- round(0.632 * n_val)
    }
    est <- tryCatch({
        fun(data)
    }, error = function(e) {
        stop("function fun/derived from fun failed with error: ", conditionMessage(e), ".\n Did you specify fun as a function or as a model object with a corresponding call and coef function?")
    })

    ## check that est is a named vector
    if (!is.vector(est) || !is.numeric(est) || is.null(names(est))) {
        stop("function fun/derived from fun did not return a named vector of coefficients")
    }
    if (parallelize) {
        requireNamespace("parallel")
        cl <- parallel::makeCluster(cores)
        parallel::clusterExport(cl, c("fun", "b", "size", "alpha", "data"), envir = environment())
        tryCatch({
            boot_est <- do.call("rbind",parallel::parLapply(cl, seq_len(b), function(i) {
                fun(data[sample(1:n_val, size, replace = FALSE), ])
            }))
        }, error = function(e) {
            stop("bootstrap (parallel) computation failed with error: ", conditionMessage(e))
        })
        parallel::stopCluster(cl)
    } else {
        pb <- utils::txtProgressBar(min = 1, max = b, style = 3, width = 50, char = "=")
        tryCatch({
            boot_est <- do.call("rbind",lapply(seq_len(b), function(i) {
                utils::setTxtProgressBar(pb, i)
                fun(data[sample(1:n_val, size, replace = FALSE), ])
            }))
        }, error = function(e) {
            stop("Bootstrap computation failed with error: ", conditionMessage(e))
        })
    }
    cat("\n")
    ## apply get_cheap_subsampling_ci for each row of boot_est, est
    tryCatch({
        res <- do.call("rbind",lapply(seq_len(length(est)), function(i) {
            get_cheap_subsampling_ci(est[i], boot_est[i, ], size, n_val, alpha)
        }))
    }, error = function(e) {
        stop("Computation of confidence intervals failed with error: ", conditionMessage(e), ". Does your function return a vector of coefficients?")
    })
    res <- cbind(Estimand = names(est),res)
    res <- list(res = res,b = b,size = size)
    if (keep_boot_estimates)
        res <- c(res,list(boot_estimates = boot_est))
    class(res) <- "cheap_subsampling"
    res
}


