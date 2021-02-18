library(shinyjs)
library(shinysky)

ui <- fluidPage(
  useShinyjs(),
  tags$body(inlineCSS(list(".shinysky-busy-indicator" = "position: absolute !important; z-index:800; "))),
  navbarPage(
    title = "Nonparametric Bounds on ATE using IVs",
    tabPanel(
      "README",
      includeMarkdown("README.md")
    ),
    tabPanel(
      "Two-Sample Bounds from Simulated Data",
      sidebarLayout(
        sidebarPanel(
          h2("Simulation Parameters"),
          wellPanel(
            numericInput(inputId = "sample_size",
                         label = "Sample Size",
                         value = 10000,
                         min = 10,
                         max = 100000),
            checkboxInput(inputId = "invalid_IVs",
                          label = "Allow Invalid IVs",
                          value = FALSE),
            checkboxInput(inputId = "include_dep_IVs",
                          label = "Include a pair of dependent variables",
                          value = FALSE)
          ),
          wellPanel(
            h3("Unmeasured Confounder"),
            # numericInput(inputId = "pU",
            #              label = HTML("P(U = 1)"),
            #              value = 0.5),
            numericInput(inputId = "U_on_X",
                         label = HTML("Effect of U on X (&gamma;<sub>U</sub>)"),
                         value = 1),
            numericInput(inputId = "U_on_Y",
                         label = HTML("Effect of U on Y (&beta;<sub>U</sub>)"),
                         value = 1)
          ),
          wellPanel(
            h3("Effect of X"),
            numericInput(inputId = "X_on_Y",
                         label = HTML("Effect of X on Y (&beta;<sub>X</sub>)"),
                         value = 2)
          ),
          wellPanel(
            h3("Instrumental Variables"),
            sliderInput(inputId = "n_indIVs", label = "Number of independent IVs to include",
                        value = 1, min = 0, max = 10, step = 1),
            radioButtons(inputId = "n_cats",
                         label = "Number of categories of Z",
                         choices = c("2", "3"), inline = TRUE),
            uiOutput(outputId = "IVs_ps"),
            h4("Effect of independent IVs on X"),
            uiOutput(outputId = "indIVs_on_X"),
            uiOutput(outputId = "indIVs_on_Y"),
            uiOutput(outputId = "depIVs")
          ),
          actionButton(inputId = "simulate_data", label = "Simulate data")
        ),
        mainPanel(
          tabsetPanel(
            id = "two_sample_bounds",
            tabPanel(
              "Introduction",
              includeMarkdown("simulation-README.md"),
              withMathJax()
            ),
            tabPanel(
              title = "Simulation Results",
              busyIndicator(),
              conditionalPanel("input.simulate_data == 0",
                               h3("Choose parameters on the left, and click 'Simulate Data' at the bottom of the sidebar panel.")),
              conditionalPanel("input.simulate_data > 0",
                               plotOutput("DAG_plotted")),
              # conditionalPanel("input.simulate_data > 0",
              #                  plotOutput("plot_results")),
              conditionalPanel(
                "input.simulate_data > 0",
                tabsetPanel(
                  tabPanel(
                    title = "Bounds Figure",
                    p(),
                    plotOutput("plot_results")
                  ),
                  tabPanel(
                    title = "Bounds",
                    p(),
                    DT::dataTableOutput("tidy_results")
                  ),
                  tabPanel(
                    title = "Simulated Data",
                    p(),
                    DT::dataTableOutput("simulated_data")
                  )
                )
              )
              #   DT::dataTableOutput("tidy_results")),
              # conditionalPanel("input.simulate_data > 0 & input.show_simulated_data",
              #                  DT::dataTableOutput("simulated_data"))
            )
          )
        )
      )
    ),
    tabPanel(
      title = "One-Sample Bounds From Two-Sample Summary Statistics",
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
                                  numericInput(inputId = "seed",
                                               label = "Specify seed",
                                               value = round(runif(1), 5)*10e4),
                                  numericInput(inputId = "dpi",
                                               label = "Figure DPI",
                                               value = 300)
                           ),
                           column(6,
                                  numericInput(inputId = "width",
                                               label = "Figure width (in)",
                                               value = 8),
                                  numericInput(inputId = "height",
                                               label = "Figure height (in)",
                                               value = 4))
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
)
