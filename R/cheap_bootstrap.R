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

get_data_and_call <- function(x, data = NULL){
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
        data <- eval(call$data, .GlobalEnv)
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
  list(data = data, n_val = n_val, x = x)
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
##' @param cores Number of cores to use for parallel computation, if parallelize = TRUE.
##' Defaults to detectCores().
##' @param data Data set to be used for the computation, if applicable. 
##' @param progress_bar Logical. Should a progress bar be displayed? Defaults to TRUE.
##' @param type Character. Type of bootstrap method. Can be either "subsampling" or "non_parametric". 
##' Defaults to "subsampling".
##' @return An object of class "cheap_bootstrap" containing 
##' the point estimates and confidence intervals.
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs <- cheap_bootstrap(x, b = 1000, data = anorexia)
##' cs
##' 
##' ## example with a model object
##' set.seed(123)
##' x <- lm(Postwt ~ Prewt + Treat + offset(Prewt), data = anorexia)
##' cs2 <- cheap_bootstrap(x, b = 1000)
##' cs2
##' 
##' ## example with a function call and parallel computation
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs3 <- cheap_bootstrap(x, b = 1000, data = anorexia, parallelize = TRUE)
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
##' cs4 <- cheap_bootstrap(ate_fit_fun, b = 100, data = dtS)
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
##' cs5 <- cheap_bootstrap(z, b = 100)
##' }
cheap_bootstrap <- function(x,
                            b = 20,
                            m_val = NULL,
                            alpha = 0.05,
                            parallelize = FALSE,
                            cores = parallel::detectCores(),
                            data = NULL,
                            progress_bar = TRUE,
                            type = "subsampling") {
  coef <- NULL
  
  data_and_call <- get_data_and_call(x, data = data)
  x <- data_and_call$x
  data <- data_and_call$data
  n_val <- data_and_call$n_val

  if (is.null(m_val) && type == "subsampling") {
    m_val <- round(0.632 * n_val)
  }
  if (type == "non_parametric"){
    m_val <- n_val
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
    parallel::clusterExport(cl, c("x", "b", "m_val", "alpha", "data", "type"), envir = environment())
    tryCatch({
      boot_est <- parallel::parSapply(cl, seq_len(b), function(i) {
        x(data[sample(1:n_val, m_val, replace = (type == "non_parametric")), ])
      })
    }, error = function(e) {
      stop("bootstrap (parallel) computation failed with error: ", conditionMessage(e))
    })
    parallel::stopCluster(cl)
  } else {
    if (progress_bar){
      pb <- utils::txtProgressBar(min = 1, max = b, style = 3, width = 50, char = "=")
    }
    tryCatch({
      boot_est <- sapply(seq_len(b), function(i) {
        if (progress_bar) utils::setTxtProgressBar(pb, i)
        x(data[sample(1:n_val, m_val, replace = (type == "non_parametric")), ])
      })
    }, error = function(e) {
      stop("bootstrap computation failed with error: ", conditionMessage(e))
    })
  }
  cat("\n")
  if (type == "non_parametric") m_val <- 1/2 * n_val ## cheap_bootstrap formula defaults to m = n/2 for non-parametric bootstrap
  ## apply get_cheap_subsampling_ci for each row of boot_est, est
  tryCatch({
    res <- sapply(seq_len(length(est)), function(i) {
      get_cheap_subsampling_ci(est[i], boot_est[i, ], m_val, n_val, alpha)
    })
  }, error = function(e) {
    stop("computation of confidence intervals failed with error: ", conditionMessage(e), ". Does your function return a vector of coefficients?")
  })
  colnames(res) <- names(est)
  res <- list(res = res, boot_estimates = boot_est, b = b, m = m_val, type = type, alpha = alpha, n = n_val)
  class(res) <- "cheap_bootstrap"
  res
}

##' Print method for cheap_bootstrap objects
##' 
##' @title Print method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Prints the point estimates and confidence intervals.
##' @export
print.cheap_bootstrap <- function(x, ...) {
  if (x$type == "subsampling") {
    cat(paste0("Cheap subsampling results for subsample size m = ", x$m, " and ", x$b, " bootstrap samples\n"))
  } else {
    cat(paste0("Cheap (non-parametric) bootstrap results for ", x$b, " bootstrap samples\n"))
  }
  print(x$res, ...)
  invisible(x)
}

## make documentation for the below function
##' Plot method for cheap_bootstrap objects
##'
##' @title Plot method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Plots the point estimates and confidence intervals as a function of the number of bootstrap samples.
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' cs <- cheap_bootstrap(x, b = 1000, data = anorexia)
##' plot(cs)
##' }
##' 
plot.cheap_bootstrap <- function(x, ...) {
  est <- as.numeric(x$res[1, ])
  if (x$type == "non_parametric") m_val <- 1/2 * x$n ## cheap_bootstrap formula defaults to m = n/2 for non-parametric bootstrap
  else m_val <- x$m
  res_b <- list()
  for (b_cur in seq_len(x$b)){
    ## apply get_cheap_subsampling_ci for each row of boot_est, est
    tryCatch({
      res <- sapply(seq_len(length(est)), function(i) {
        get_cheap_subsampling_ci(est[i], x$boot_est[i, seq_len(b_cur)], m_val, x$n, x$alpha)
      })
    }, error = function(e) {
      stop("computation of confidence intervals failed with error: ", conditionMessage(e), ". Does your function return a vector of coefficients?")
    })
    colnames(res) <- colnames(x$res)
    ## widen res, so that the columns 
    res <- data.frame(lapply(data.frame(t(res)), unlist), stringsAsFactors=FALSE)
    res$type <- colnames(x$res)
    res$b <- b_cur
    res_b[[b_cur]] <- res
  }
  res_b <- do.call(rbind, res_b)
  ggplot2::ggplot(data = res_b, ggplot2::aes(x = b, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.2, ggplot2::aes(ymin = lower_b, ymax = upper_b)) +
    ggplot2::facet_wrap(~type, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Number of bootstrap samples") +
    ggplot2::theme(legend.position = "bottom")
}

get_shiny_panel <- function(n_val){
  shiny::sidebarPanel(
    shiny::sliderInput("b", "Number of bootstrap samples", min = 1, max = 1000, value = 20, step = 1),
    shiny::sliderInput("m", "Subsample size", min = round(0.2* n_val), max = n_val, value = round(0.632 * n_val), step = 1),
    shiny::sliderInput("alpha", "Significance level", min = 0.01, max = 0.1, value = 0.05),
    shiny::sliderInput("cores", "Number of cores", min = 1, max = parallel::detectCores(), value = parallel::detectCores() - 1, step = 1),
    shiny::checkboxInput("parallelize", "Parallelize computation", value = FALSE),
    shiny::selectInput("type", "Bootstrap type", choices = c("subsampling", "non_parametric"), selected = "subsampling"),
    shiny::actionButton("run", "Run")
  )
}

##' Shiny app for cheap subsampling
##' 
##' This function creates a shiny app for the cheap subsampling method.
##' The user can vary the number of bootstrap samples and the subsample size.
##' 
##' @title Shiny app for cheap subsampling
##' @param x A function that returns a vector of coefficients 
##' or a model object which saves the call. In the second case, a
##' coef needs to be defined for 'class(x)'.
##' @param data Data set to be used for the computation, if applicable. 
##' @export
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' shiny_cheap_bootstrap(x, data = anorexia)
shiny_cheap_bootstrap <- function(x, data = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for this function")
  }
  
  data_and_call <- get_data_and_call(x, data = data)
  x <- data_and_call$x
  data <- data_and_call$data
  n_val <- data_and_call$n_val
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("Cheap bootstrap"),
    shiny::sidebarLayout(
      get_shiny_panel(n_val),
      shiny::mainPanel(
        shiny::tableOutput("table")
      )
    )
  )
  server <- function(input, output) {
    ## on event run, run the cheap_bootstrap function
    shiny::observeEvent(input$run, {
      output$table <- shiny::renderTable({
        cheap_bootstrap(
          x,
          b = input$b,
          m_val = input$m,
          alpha = input$alpha,
          parallelize = input$parallelize,
          cores = input$cores,
          data = data,
          type = input$type,
          progress_bar = FALSE
        )$res
      })
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}

## make a shiny app for the plot method
##' Shiny app for cheap subsampling plot
##'
##' @title Shiny app for cheap subsampling
##' @param x A function that returns a vector of coefficients 
##' or a model object which saves the call. In the second case, a
##' coef needs to be defined for 'class(x)'.
##' @param data Data set to be used for the computation, if applicable. 
##' @export
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' x <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' shiny_cheap_bootstrap_plot(x, data = anorexia)
##' 
shiny_cheap_bootstrap_plot <- function(x, data = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for this function")
  }
  
  data_and_call <- get_data_and_call(x, data = data)
  x <- data_and_call$x
  data <- data_and_call$data
  n_val <- data_and_call$n_val
  
  ui <- shiny::fluidPage(
      get_shiny_panel(n_val),
      shiny::plotOutput("plot")
    )
  server <- function(input, output) {
    ## on event run, run the cheap_bootstrap function
    shiny::observeEvent(input$run, {
      output$plot <- shiny::renderPlot(
        plot.cheap_bootstrap(cheap_bootstrap(
          x,
          b = input$b,
          m_val = input$m,
          alpha = input$alpha,
          parallelize = input$parallelize,
          cores = input$cores,
          data = data,
          type = input$type,
          progress_bar = FALSE
        ))
      )
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}
