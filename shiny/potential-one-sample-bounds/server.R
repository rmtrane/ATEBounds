library(shiny)
library(ATEBounds)
# library(dagitty)
# library(ggdag)
# library(ggrepel)
library(tidyverse)

if(FALSE){
  input <- list(theta0 = 0.2, theta1 = 0.2, theta2 = 0.2,
                gamma0 = 0.2, gamma1 = 0.2, gamma2 = 0.2,
                x_mono = TRUE, seed = 26971,
                n_joints = 50)
}


server <- function(input, output){
  ## One-Sample Bounds From Two-Sample Summary Statistics
  observeEvent(input$sim, {
    thetas <- runif(3)

    gammas <- vector(length = length(thetas))


    if(input$x_mono == "No"){
      gammas[1] <- runif(1)

      gammas[2] <- runif(1,
                         min = max(0,
                                   gammas[1] - thetas[1] - thetas[2], # row 3
                                   gammas[1] + thetas[1] + thetas[2] - 2), # row 5
                         max = min(1,
                                   gammas[1] - thetas[1] - thetas[2] + 2, # row 11
                                   gammas[1] + thetas[1] + thetas[2])) # row 16

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
    } else {
      # ATEBounds:::matrices_from_polymake %>% filter(x_monotone, n_z_levels == 3, data_format == "bivariate", !y_monotone) %>% pull(matrix) %>% .[[1]] %>% filter(alpha == 0) %>% select(where(~sum(abs(.x)) > 0)) %>% filter(rowSums(abs(.)) > 1)

      thetas <- sort(thetas)

      gammas[1] <- runif(1)

      gammas[2] <- runif(1,
                         min = max(0,
                                   gammas[1] + thetas[1] - thetas[2]), # row 8
                         max = min(1,
                                   gammas[1] - thetas[1] + thetas[2] # row 4
                                   ))

      gammas[3] <- runif(1,
                         min = max(0,
                                   gammas[2] + thetas[2] - thetas[3], # row 3
                                   gammas[2] - gammas[1]), # row 5
                         max = min(1,
                                   gammas[2] - thetas[2] + thetas[3], # row 2
                                   1 - gammas[1] + gammas[2] # row 7
                                   ))
    }

    thetas <- round(thetas, 4)
    gammas <- round(gammas, 4)

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
    updateTabsetPanel(inputId = "main_panels", selected = "one_sample_bounds")
  })

  observeEvent({
    input$run
    input$main_panels
  }, {
    if(input$main_panels == "one_sample_bounds" && input$run > 0){

    set.seed(input$seed)
    cat(input$seed, "\n")

    RVs$two_sample_bounds <- get_bounds(gammas = c(input$gamma0, input$gamma1, input$gamma2),
                                        thetas = c(input$theta0, input$theta1, input$theta2),
                                        stop = FALSE, warning = TRUE,
                                        x_mono = input$x_mono == "Yes")

    if(!RVs$two_sample_bounds$constraints_violated){

      RVs$pot_covs <- potential_covs(thetas = c(input$theta0, input$theta1, input$theta2),
                                     gammas = c(input$gamma0, input$gamma1, input$gamma2),
                                     x_mono = input$x_mono == "Yes")

      RVs$joint_samples <- sample_joint_probs(pot_covs = RVs$pot_covs,
                                              n = input$n_joints,
                                              max_rejections = ceiling(input$n_joints/10),
                                              x_mono = input$x_mono == "Yes",
                                              return_bounds = TRUE) %>%
        unnest(joint) %>%
        unnest_wider(bounds) %>%
        arrange(lower) %>%
        mutate(contains_zero = if_else(lower > 0 | upper < 0, "Does Not Overlap Zero", "Overlaps Zero"),
               id = row_number() / input$n_joints)

      RVs$pretty_plot <- ggplot(filter(RVs$joint_samples, !is.na(contains_zero)),
                                aes(xmin = lower, xmax = upper,
                                    y = id,
                                    color = contains_zero,
                                    alpha = contains_zero)) +
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
        geom_vline(xintercept = 0) +
        geom_errorbar(aes(linetype = contains_zero)) +
        geom_vline(data = data.frame(),
                   aes(xintercept = RVs$two_sample_bounds$interval,
                       linetype = rep("Two-Sample Bounds", 2),
                       color = rep("Two-Sample Bounds", 2))
        ) +
        scale_x_continuous(limits = c(-1, 1),
                           breaks = c(-1, -0.5, 0, 0.5, 1)) +
        scale_color_manual(values = c("black", "blue", "red"),
                           breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does Not Overlap Zero")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid"),
                              breaks = c("Two-Sample Bounds", "Overlaps Zero", "Does Not Overlap Zero")) +
        scale_alpha_manual(values = c(0.5, 1),
                           breaks = c("Overlaps Zero", "Does Not Overlap Zero")) +
        guides(
          alpha = "none"
        ) +
        labs(
          x = "ATE",
          y = "",
          color = "",
          fill = "",
          linetype = "",
          caption = paste0("Proportion of one-sample bounds not overlapping 0: ", 100*round(mean(RVs$joint_samples$contains_zero == "Does Not Overlap Zero", na.rm = TRUE), digits = 3), "%.")
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
                 dpi = input$dpi,
                 width = input$width, height = input$height)
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
            plotOutput("prettyPlot", height = input$plot_height)
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
    }
  })
}
