get_shiny_panel <- function(n_val){
  shiny::sidebarPanel(
    shiny::sliderInput("b", "Number of bootstrap samples", min = 1, max = 1000, value = 20, step = 1),
    shiny::sliderInput("size", "Subsample size", min = round(0.2* n_val), max = n_val, value = round(0.632 * n_val), step = 1),
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
##' @param fun A function that returns a vector of coefficients 
##' or a model object which saves the call. In the second case, a
##' coef needs to be defined for 'class(fun)'.
##' @param data Data set to be used for the computation, if applicable. 
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' fun <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' shiny_cheap_bootstrap(fun, data = anorexia)
##' }
shiny_cheap_bootstrap <- function(fun, data = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for this function")
  }
  
  data_and_call <- get_data_and_call(fun, data = data)
  fun <- data_and_call$fun
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
          fun,
          b = input$b,
          size = input$size,
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
##' @param fun A function that returns a vector of coefficients 
##' or a model object which saves the call. In the second case, a
##' coef needs to be defined for 'class(fun)'.
##' @param data Data set to be used for the computation, if applicable. 
##' @export
##' @examples
##' \dontrun{
##' utils::data(anorexia, package = "MASS")
##' ## example with a function call
##' set.seed(123)
##' fun <- function(d) coef(lm(Postwt ~ Prewt + Treat + offset(Prewt), data = d))
##' shiny_cheap_bootstrap_plot(fun, data = anorexia)
##' }
##' 
shiny_cheap_bootstrap_plot <- function(fun, data = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for this function")
  }
  
  data_and_call <- get_data_and_call(fun, data = data)
  fun <- data_and_call$fun
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
          fun,
          b = input$b,
          size = input$size,
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
