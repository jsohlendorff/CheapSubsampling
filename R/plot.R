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
##' cs <- cheap_bootstrap(x, b = 10, data = anorexia)
##' plot(cs)
##' }
##' 
plot.cheap_bootstrap <- function(x, ...) {
  b <- estimate <- lower_b <- upper_b <- NULL 
  est <- as.numeric(x$res[1, ])
  if (x$type == "non_parametric") size <- 1/2 * x$n ## cheap_bootstrap formula defaults to m = n/2 for non-parametric bootstrap
  else size <- x$size
  res_b <- list()
  for (b_cur in seq_len(x$b)){
    ## apply get_cheap_subsampling_ci for each row of boot_est, est
    tryCatch({
      res <- sapply(seq_len(length(est)), function(i) {
        get_cheap_subsampling_ci(est[i], x$boot_est[i, seq_len(b_cur)], size, x$n, x$alpha)
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
