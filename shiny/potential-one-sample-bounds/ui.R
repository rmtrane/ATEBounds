library(shinyjs)
library(shinysky)

ui <- fluidPage(
  useShinyjs(),
  tags$body(inlineCSS(list(".shinysky-busy-indicator" = "position: absolute !important; z-index:800; "))),
  titlePanel("One-Sample Bounds From Two-Sample Summary Statistics"),
  sidebarLayout(
    sidebarPanel(
      column(12,
             radioButtons(inputId = "x_mono",
                          label = "Assume P(X = 1 | Z = z, U) < P(X = 1 | Z = z+1, U)?",
                          choices = list("Yes",
                                         "No"),
                          selected = "No"),
             sliderInput(inputId = "n_joints",
                         label = "Sampling Iterations",
                         min = 1, max = 1000,
                         step = 1,
                         value = 50)
      ),
      fluidRow(
        column(4,
               numericInput("theta0",
                            "P(X = 1 | Z = 0)",
                            min = 0,
                            max = 1,
                            value = 0.2)),
        column(4,
               numericInput("theta1",
                            "P(X = 1 | Z = 1)",
                            min = 0,
                            max = 1,
                            value = 0.2)),
        column(4,
               numericInput("theta2",
                            "P(X = 1 | Z = 2)",
                            min = 0,
                            max = 1,
                            value = 0.2))
      ),
      fluidRow(
        column(4,
               numericInput("gamma0",
                            "P(Y = 1 | Z = 0)",
                            min = 0,
                            max = 1,
                            value = 0.2)),
        column(4,
               numericInput("gamma1",
                            "P(Y = 1 | Z = 1)",
                            min = 0,
                            max = 1,
                            value = 0.2)),
        column(4,
               numericInput("gamma2",
                            "P(Y = 1 | Z = 2)",
                            min = 0,
                            max = 1,
                            value = 0.2)
        ),
        column(12,
               actionButton(inputId = "sim",
                            label = "Create Random Values"),
               actionButton(inputId = "run",
                            label = "GO!")
        )
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = "main_panels",
        tabPanel(
          "Introduction",
          includeMarkdown("one-sample-README.md"),
          withMathJax()
        ),
        tabPanel(
          value = "one_sample_bounds",
          "Potential One-Sample Bounds",
          busyIndicator(),
          conditionalPanel(
            "input.run == 0",
            h2("Specify probabilities on the left and click the 'GO!' button")
          ),
          conditionalPanel(
            "input.run != 0",
            uiOutput("error_or_plot")
          ),
          fluidRow(
            column(12,
                   checkboxInput(inputId = "advanced",
                                 label = "Advanced"),
                   conditionalPanel(
                     condition = "input.advanced",
                     fluidRow(
                       column(6,
                              h2("In app options:"),
                              numericInput(inputId = "seed",
                                           label = "Specify seed",
                                           value = round(runif(1), 5)*10e4),
                              sliderInput("plot_height",
                                          label = "Plot height (pixels)",
                                          min = 400,
                                          max = 1200,
                                          value = 400)
                       ),
                       column(6,
                              h2("Download options:"),
                              numericInput(inputId = "width",
                                           label = "Figure width (in)",
                                           value = 8),
                              numericInput(inputId = "height",
                                           label = "Figure height (in)",
                                           value = 4),
                              numericInput(inputId = "dpi",
                                           label = "Figure DPI",
                                           value = 300)
                       )
                     )
                   )
            )
          )
        )
      )
    )
  )
)
