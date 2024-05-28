get_data_and_call <- function(fun, data = NULL){
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
        data <- eval(call$data, .GlobalEnv)
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
    stop("Data needs to be a data frame")
  }
  n_val <- nrow(data)
  list(data = data, n_val = n_val, fun = fun)
}
