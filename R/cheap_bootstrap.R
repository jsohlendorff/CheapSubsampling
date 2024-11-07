## retrieve cheap subsampling confidence interval
get_conf_int <- function(est, est_boot, size, n_dat, alpha, type = "subsampling") {
  b_max <- length(est_boot)
  res <- list()
  for (b in seq_len(b_max)){
    b_est <- est_boot[seq_len(b_max)]
    if (type == "subsampling") {
      sub_factor <- sqrt((size) / (n_dat - size))
    } else {
      sub_factor <- 1
    }
    s_val <- sqrt(mean((est - b_est)^2))
    tq <- stats::qt(1 - alpha / 2, df = b)
    res[[b]] <- data.frame(
      estimate = est,
      cheap_lower = est - tq * sub_factor * s_val,
      cheap_upper = est + tq * sub_factor * s_val,
      b = b
    )
  }
  do.call(rbind, res)
}

##' Method implementing the cheap subsampling method for confidence intervals
##'
##' Given a model object or a function that returns a vector of coefficients
##' and a data set,
##' this function computes confidence intervals using
##' the cheap subsampling method.
##'
##' @title Cheap subsampling
##' @param fun A function that returns a data.frame with the point estimates given
##' in est_col_name and the corresponding parameter names in par_col_names.
##' @param data Data set to be used for the computation.
##' @param par_col_names Character. Column names of the parameter names.
##' @param est_col_name Character. Column name of the point estimates.
##' @param b Number of bootstrap samples.
##' @param size Subsample size. Defaults to 0.632 * nrow(data).
##' @param alpha Significance level. Defaults to 0.05.
##' @param type Character. Type of bootstrap method.
##' Can be either "subsampling" or "non_parametric".
##' Defaults to "subsampling".
##' @param parallelize Logical. Should the computation be parallelized?
##' @param cores Number of cores to be used for parallel computation.
##' @param keep_estimates Logical. Should the bootstrapped estimates be kept?
##' @return An object of class "cheap_bootstrap" containing
##' the point estimates and confidence intervals.
##' @export
##' @examples
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' library(broom)
##' library(dplyr)
##' fun <- function(data) {
##'   lm(formula = Postwt ~ Prewt + Treat + offset(Prewt), data = data) %>%
##'   tidy() 
##' }
##' cs <- cheap_bootstrap(fun = fun, 
##'                       b = 20, 
##'                       data = anorexia,
##'                       est_col_name = "estimate",
##'                       par_col_names = "term")
##' cs
##' summary(cs)
##' plot(cs)
##'
##' ## example with a function call and parallel computation
##' \dontrun{
##' set.seed(123)
##' ## note the function needs to load the packages needed
##' ## for the computation if parallelized
##' cs <- cheap_bootstrap(fun = fun, 
##'                       b = 20, 
##'                       data = anorexia, 
##'                       est_col_name = "estimate",
##'                       par_col_names = "term", 
##'                       parallelize = TRUE, 
##'                       cores = 2)
##' }
cheap_bootstrap <- function(fun,
                            data,
                            size = NULL,
                            type = "subsampling",
                            b = 20,
                            est_col_name = "estimate",
                            par_col_names = "parameter",
                            alpha = 0.05,
                            parallelize = FALSE,
                            cores = 1,
                            keep_estimates = TRUE) {
  estimate <- b_est <- ..keep <- ..keep_names <- NULL
  ## Check if the data is a data.frame
  if (!inherits(data, "data.frame")) {
    stop("data must inherit data.frame")
  }
  n_dat <- nrow(data)
  
  ## If parallelize, check that namespace is loaded
  if (parallelize){
    requireNamespace("parallel")
  }

  ## Check that all relevant arguments are of length 1
  arg_check <- c(
    "b",
    "alpha",
    "type",
    "parallelize",
    "cores",
    "keep_estimates"
  )
  for (var_name in arg_check) {
    if (length(get(eval(var_name))) != 1) {
      stop(paste(var_name, "must be of length 1"))
    }
  }

  if (is.null(size) && type == "subsampling") {
    size <- round(0.632 * n_dat)
  } else if (type == "non_parametric") {
    size <- n_dat
  }

  est <- tryCatch(
    {
      fun(data)
    },
    error = function(e) {
      stop(
        "function fun/derived from fun failed with error: ",
        conditionMessage(e),
        ".\n Did you specify fun as a function or as a model object with a corresponding call and coef function?"
      )
    }
  )
  ## Check that est inherits data.frame
  if (!inherits(est, "data.frame")) {
    stop("function fun/derived from fun must return a data.frame")
  }
  
  est <- data.table::as.data.table(est)
  
  ## Check that est_col_name and par_col_names are in the names of columns of est
  if (!(est_col_name %in% colnames(est)) || 
      !(all(par_col_names %in% colnames(est)))) {
    stop("est_col_name and par_col_names must be in the names of columns of est")
  }
  
  ## Check that that the column specified by est_col_name is numeric
  if (!all(sapply(est[[est_col_name]], is.numeric))) {
    stop("The column specified by est_col_name must be numeric")
  }
  
  ## Check that the column specified by par_col_names is character or factor
  if (sum(est[, sapply(.SD, function(x) !is.character(x) && !is.factor(x)), .SDcols = par_col_names]) > 0) {
    stop("The column specified by par_col_names must be character or factor")
  }
  keep_names <- c(par_col_names, est_col_name)
  est <- est[, ..keep_names]
  
  boot_fun <- function(i) {
    set.seed(seeds[i])
    tryCatch(
      {
        x<- fun(data[
          sample(1:n_dat,
                 size,
                 replace = (type == "non_parametric")
          ), ,
          drop = FALSE
        ])
        x$b <- i
        x
      },
      error = function(e) {
        stop(
          "Bootstrap computation failed with error: ",
          conditionMessage(e),
          " for iteration ", i," with the seed ", seeds[i]
        )
      }
    )
  }
  
  ## sample b seeds
  seeds <- sample.int(1e+09, b)
  
  if (parallelize) {
    results <- pbmcapply::pbmclapply(X = seq_len(b),
                                     FUN = boot_fun,
                                     mc.cores = cores)
  } else {
    results <- lapply(seq_len(b), boot_fun)
  }
  tryCatch(
    {
      keep <- c(par_col_names, est_col_name, "b")
      boot_est <- data.table::as.data.table(do.call("rbind", results))[, ..keep]
    },
    error = function(e) {
      stop("Could not bind the results. Did you return a vector of coefficients?")
    }
  )
  
  ## Rename the column by est_col_name to b_est
  ## and rename the column from est to estimate
  data.table::setnames(boot_est, est_col_name, "b_est")
  data.table::setnames(est, est_col_name, "estimate")
  
  ## Merge boot_est with by par_col_names
  boot_est <- merge(est, boot_est, by = par_col_names, all.x = TRUE)
  boot_est[, c("size", "n_dat", "alpha", "type") := list(size, n_dat, alpha, type)]

  tryCatch(
    {
      ## Apply get_conf_int for each column specified by par_col_names
      res <- boot_est[, get_conf_int(estimate[1], b_est, size[1], n_dat[1], alpha[1], type[1]), by = par_col_names]
    },
    error = function(e) {
      stop(
        "Computation of confidence intervals failed with error: ",
        conditionMessage(e),
        ". Does your function return a vector of coefficients?"
      )
    }
  )
  
  res <- list(
    res = res,
    b = b,
    size = size,
    type = type,
    alpha = alpha,
    n = n_dat,
    par_col_names = par_col_names
  )
  if (keep_estimates) {
    res <- c(res, list(boot_estimates = boot_est))
  }
  class(res) <- "cheap_bootstrap"
  res
}

##' Summary method for cheap_subsampling objects
##'
##' @title Summary method for cheap_subsampling objects
##' @param object An object of class "cheap_subsampling"
##' @param print Logical. Should the summary be printed? Defaults to TRUE.
##' @param ... Not applicable.
##' @return Summary table with the point estimates and confidence intervals.
##' @export
summary.cheap_bootstrap <- function(object, print = TRUE, ...) {
  if (print[[1]] == TRUE) {
    print(object, ...)
  }
  invisible(object$res)
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
    cat(
      paste0(
        "Cheap subsampling results for subsample size m = ",
        x$size,
        " and ",
        x$b,
        " bootstrap samples\n"
      )
    )
  } else {
    cat(
      paste0(
        "Cheap (non-parametric) bootstrap results for ",
        x$b,
        " bootstrap samples\n"
      )
    )
  }
  print(x$res, ...)
  invisible(x)
}

##' Plot method for cheap_bootstrap objects
##'
##' @title Plot method for cheap_bootstrap objects
##' @param x An object of class "cheap_bootstrap"
##' @param ... Not applicable.
##' Plots the point estimates and confidence intervals
##' as a function of the number of bootstrap samples.
##' @export
plot.cheap_bootstrap <- function(x, ...) {
  form_par_wrap <- stats::as.formula(paste0("~", paste0(x$par_col_names, collapse = "+")))
  b <- estimate <- cheap_lower <- cheap_upper <- NULL
  ggplot2::ggplot(data = x$res, 
                  ggplot2::aes(x = b, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      alpha = 0.2,
      ggplot2::aes(
        ymin = cheap_lower,
        ymax = cheap_upper
      )
    ) +
    ggplot2::facet_wrap(form_par_wrap, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Number of bootstrap samples") +
    ggplot2::ylab(paste0(
      "Cheap ",
      ifelse(x$type == "subsampling", "subsampling", "bootstrap"),
      " confidence intervals"
    ))
}