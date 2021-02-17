#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ATEBounds)
library(shinyjs)
library(shinysky)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Setup busy indicator
    useShinyjs(),
    tags$body(inlineCSS(list(".shinysky-busy-indicator" = "position: absolute !important; z-index:800; "))),
    navbarPage(
        # Application title
        "One-Sample Bounds From Two-Sample Summary Statistics",

        tabPanel(
            "Potential One-Sample Bounds",
            # Sidebar with a slider input for number of bins
            sidebarLayout(
                sidebarPanel(
                    column(12,
                           numericInput(inputId = "seed",
                                        label = "Seed to use",
                                        value = round(runif(1), 5)*10e4),
                           radioButtons(inputId = "x_mono",
                                        label = "Assume P(X = 1 | Z = z) < P(X = 1 | Z = z+1)?",
                                        choices = list("Yes" = TRUE,
                                                       "No" = FALSE),
                                        selected = FALSE),
                           sliderInput(inputId = "n_joints",
                                       label = "Number of joint distributions to sample",
                                       min = 1, max = 1000,
                                       step = 1,
                                       value = 100)
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
                        )
                    ),
                    fluidRow(
                        column(
                            6,
                            actionButton(inputId = "run",
                                         label = "GO!")
                        ),
                        column(
                            6,
                            actionButton(inputId = "sim",
                                         label = "Create Random Values")
                        )
                    ),
                    textOutput("go")
                ),

                # Show a plot of the generated distribution
                mainPanel(
                    busyIndicator(),
                    uiOutput("error_or_plot")
                )
            )
        ),
        tabPanel(
            "README",
            includeMarkdown("README.md"),
            withMathJax()
        )
    )
)


if(FALSE){
    input <- list(theta0 = 0.2, theta1 = 0.2, theta2 = 0.2,
                  gamma0 = 0.2, gamma1 = 0.2, gamma2 = 0.2,
                  x_mono = TRUE, seed = 26971,
                  n_joints = 100)
}

# Define server logic required to draw a histogram
server <- function(input, output) {

    RVs <- list(
        pot_covs = NULL,
        joint_samples = NULL,
        two_sample_bounds = NULL,
        continue = 0
    )

    observeEvent(input$sim, {
        thetas <- runif(3)

        gammas <- vector(length = length(thetas))


        gammas[1] <- runif(1)

        gammas[2] <- runif(1,
                           min = max(0,
                                     gammas[1] - thetas[1] - thetas[2], # row 3
                                     gammas[1] + thetas[1] + thetas[2] - 2), # row 5
                           max = min(1,
                                     gammas[1] - thetas[1] - thetas[2] + 2, # row 11
                                     gammas[1] + thetas[1] + thetas[2]) # row 16
        )

        gammas[3] <- runif(1,
                           min = max(0,
                                     gammas[1] - thetas[1] - thetas[3], # row 4
                                     gammas[1] + thetas[1] + thetas[3] - 2, # row 6
                                     gammas[2] - thetas[2] - thetas[3], # row 12
                                     gammas[2] + thetas[2] + thetas[3] - 2), # row 9
                           max = min(1,
                                     gammas[2] - thetas[2] - thetas[3] + 2, # row 10
                                     gammas[2] + thetas[2] + thetas[3], # row 13
                                     gammas[1] - thetas[1] - thetas[3] + 2, # row 17
                                     gammas[1] + thetas[1] + thetas[3])) # row 15

        updateNumericInput(inputId = "theta0", value = thetas[1])
        updateNumericInput(inputId = "theta1", value = thetas[2])
        updateNumericInput(inputId = "theta2", value = thetas[3])
        updateNumericInput(inputId = "gamma0", value = gammas[1])
        updateNumericInput(inputId = "gamma1", value = gammas[2])
        updateNumericInput(inputId = "gamma2", value = gammas[3])
    })

    observeEvent(input$n_joints, {
        if(input$n_joints > 200){
            showModal(
                modalDialog(
                    title = "Large n",
                    p("Number of samples requested is large. This might take a while."),
                    footer = tagList(
                        actionButton("reset_n_joints",
                                     label = "Cancel"),
                        modalButton("Continue")
                        # actionButton("continue",
                        #              label = "Continue")
                    )
                )
            )
        }
    })

    observeEvent(input$reset_n_joints, {
        removeModal()
        updateSliderInput(inputId = "n_joints", value = 200)
    })

    observeEvent(input$continue, {
        removeModal()
    })

    observeEvent(input$run, {
        set.seed(input$seed)
        cat(input$seed, "\n")

        RVs$two_sample_bounds <- get_bounds(gammas = c(input$gamma0, input$gamma1, input$gamma2),
                                            thetas = c(input$theta0, input$theta1, input$theta2),
                                            stop = FALSE, warning = TRUE,
                                            x_mono = input$x_mono)

        if(!RVs$two_sample_bounds$constraints_violated){

            RVs$pot_covs <- potential_covs(thetas = c(input$theta0, input$theta1, input$theta2),
                                           gammas = c(input$gamma0, input$gamma1, input$gamma2),
                                           x_mono = input$x_mono)

            RVs$joint_samples <- sample_joint_probs(pot_covs = RVs$pot_covs,
                                                    n = input$n_joints,
                                                    max_rejections = ceiling(input$n_joints/10),
                                                    x_mono = input$x_mono,
                                                    return_bounds = TRUE) %>%
                unnest(joint) %>%
                unnest_wider(bounds) %>%
                arrange(lower) %>%
                mutate(contains_zero = if_else(lower > 0 | upper < 0, "Does Not Overlap Zero", "Overlaps Zero"),
                       id = row_number() / input$n_joints)

            RVs$pretty_plot <- ggplot(RVs$joint_samples,
                                      aes(xmin = lower, xmax = upper,
                                          y = id,
                                          color = contains_zero)) +
                geom_rect(data = data.frame(x_min = rep(0.2, 3),
                                            x_max = rep(0.3, 3),
                                            y_min = rep(0.1, 3),
                                            y_max = rep(0.2, 3),
                                            contains_zero = c("Two-Sample Bounds", "Overlaps Zero", "Does Not Overlap Zero")),
                          inherit.aes = FALSE,
                          aes(xmin = x_min,
                              xmax = x_max,
                              ymin = y_min,
                              ymax = y_max,
                              linetype = contains_zero,
                              color = contains_zero),
                          alpha = 0,
                          size = 0) +
                geom_errorbar(aes(linetype = contains_zero)) +
                geom_vline(data = data.frame(),
                           aes(xintercept = RVs$two_sample_bounds$interval,
                               linetype = rep("Two-Sample Bounds", 2),
                               color = rep("Two-Sample Bounds", 2))
                ) +
                scale_x_continuous(limits = c(-1, 1),
                                   breaks = c(-1, -0.5, 0, 0.5, 1)) +
                scale_color_manual(values = c("black", "grey50", "red"),
                                   breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does Not Overlap Zero")) +
                scale_linetype_manual(values = c("dashed", "solid", "solid"),
                                      breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does Not Overlap Zero")) +
                labs(
                    x = "ATE",
                    y = "",
                    color = "",
                    fill = "",
                    linetype = "",
                    caption = paste0("Proportion of one-sample bounds not overlapping 0: ", 100*round(mean(RVs$joint_samples$contains_zero != "Overlaps Zero"), digits = 3), "%.")
                ) +
                theme_bw() +
                theme(axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      text = element_text(size = 20),
                      legend.position = "top")

            output$prettyPlot <- renderPlot({
                print(RVs$pretty_plot)
            })

            output$rds_data <- downloadHandler(
                filename = function() {
                    paste0("one-sample-bounds-", Sys.Date(), ".Rds")
                },
                content = function(file) {
                    write_rds(x = RVs$joint_samples %>%
                                  mutate(two_sample_lower = RVs$two_sample_bounds$interval[["lower"]],
                                         two_sample_upper = RVs$two_sample_bounds$interval[["upper"]]),
                              file = file)
                }
            )

            output$pretty_plot <- downloadHandler(
                filename = function() {
                    paste0("one-sample-bounds-", Sys.Date(), ".png")
                },
                content = function(file) {
                    ggsave(filename = file,
                           RVs$pretty_plot,
                           dpi = 300,
                           width = 8, height = 4)
                }
            )
        }

        output$error_or_plot <- renderUI({
            if(RVs$two_sample_bounds$constraints_violated){
                h2("The specified values of P(X = 1 | Z = z) and P(Y = 1 | Z = z) violate one or more of the two-sample IV constraints.")
            } else {
                fluidRow(
                    column(
                        12,
                        plotOutput("prettyPlot")
                    ),
                    column(
                        6,
                        downloadButton("rds_data", label = "Save one-sample data")
                    ),
                    column(
                        6,
                        downloadButton("pretty_plot", label = "Save plot")
                    )
                )
            }
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
