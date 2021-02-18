library(shiny)
library(ATEBounds)
library(dagitty)
library(ggdag)
library(ggrepel)
library(tidyverse)

# source("../R/misc.R")

plot_bounds <- function(results, sim_data_coefficients){

  # tidied_output <- my_bds %>%
  #   mutate(tidy_output = map(bpbounds, tidy_bpbounds),
  #          `n Z categories` = map_dbl(bpbounds, function(x) x$nzcats)) %>%
  #   unnest(tidy_output)

  for_plot <- results %>%
    filter(!constraints_violated) %>%
    mutate(IV = paste0("Z", p)) %>%
    left_join(
      sim_data_coefficients %>%
        filter(str_detect(effect, "Z[0-9]+_on_X")) %>%
        separate(effect, into = c("IV", "X"), sep = "_on_") %>%
        select(-X)
    ) %>%
    arrange(p)

  IVs <- for_plot$IV
  IV_coefs <- for_plot$coef

  true_X_on_Y <- sim_data_coefficients %>% filter(effect == "X_on_Y") %>% pull(coef)
  true_Y_int <- sim_data_coefficients %>% filter(effect == "Yintercept") %>% pull(coef)
  #true_U_on_Y <- sim_data_coefficients %>% filter(effect == "U_on_Y") %>% pull(coef)

  true_effect <- plogis(true_Y_int + true_X_on_Y) - plogis(true_Y_int)

  ggplot(for_plot,
         aes(xmin = lower, xmax = upper, y = (p - 0.5) / length(IVs))) +
    geom_vline(xintercept = true_effect, color = "red", linetype = "dashed") +
    geom_errorbar(width = 0.1) +
    scale_y_continuous(name = "IV",
                       limits = c(0, 1),
                       labels = IVs,
                       breaks = (seq_along(IVs) - 0.5) / length(IVs),
                       sec.axis = sec_axis(trans = function(x) x,
                                           breaks = (seq_along(IVs) - 0.5) / length(IVs),
                                           labels = IV_coefs,
                                           name = "Effect on X")) +
    scale_x_continuous(limits = c(-1,1)) +
    theme_bw()
}

if(FALSE){
  input <- list(theta0 = 0.2, theta1 = 0.2, theta2 = 0.2,
                gamma0 = 0.2, gamma1 = 0.2, gamma2 = 0.2,
                x_mono = TRUE, seed = 26971,
                n_joints = 50)
}


server <- function(input, output){
  output$indIVs_on_X <- renderUI({
    if(input$n_indIVs > 0){
      lapply(1:input$n_indIVs,
             function(i){
               numericInput(inputId = paste0("Z", i, "_on_X"),
                            label = HTML(paste0("Z<sub>", i, "</sub> (&gamma;<sub>", i, "</sub>)")),
                            value = 0.2*i)
             })
    }
  })

  output$indIVs_on_Y <- renderUI({
    if(input$n_indIVs > 0 & input$invalid_IVs){
      lapply(0:input$n_indIVs,
             function(i){
               if(i < 1){
                 h3("Effect of independent IVs on Y")
               } else {
                 numericInput(inputId = paste0("Z", i, "_on_Y"),
                              label = HTML(paste0("Z<sub>", i, "</sub> (&beta;<sub>", i, "</sub>)")),
                              value = 0)
               }
             })
    }
  })


  # output$IVs_ps0 <- renderUI({
  #   map2(1:(as.numeric(input$n_indIVs)),
  #        1:(as.numeric(input$n_cats)-1),
  #        function(x = .x, y = .y){
  #
  #        })
  #   lapply(1:(as.numeric(input$n_cats)-1),
  #          function(i){
  #            numericInput(inputId = paste0("IVs_ps", i, "_", k),
  #                         label = HTML(glue("P(Z<sub>{k}</sub> = {i-1})")),
  #                         value = round(1/as.numeric(input$n_cats), digits = 3),
  #                         min = 0,
  #                         max = 1)
  #          })
  # })

  output$IVs_ps <- renderUI({
    expand_grid(x = 1:as.numeric(input$n_indIVs),
                y = 1:(as.numeric(input$n_cats)-1)) %>%
      mutate(outs = map2(x, y,
                         .f = function(x,y){
                           numericInput(inputId = paste0("IVs_ps", x, "_", y), #P(Z_x = y - 1)
                                        label = HTML(paste0("P(Z<sub>", x, "</sub> = ", y-1, ")")),
                                        value = ifelse(input$n_cats == 2, 0.5,
                                                       c(0.25, 0.5, 0.25)[y]), #round(1/as.numeric(input$n_cats), digits = 3),
                                        min = 0,
                                        max = 1)
                         })
      ) %>% pull(outs) %>% return()
  })

  output$depIVs <- renderUI({
    if(input$include_dep_IVs){
      outs <- list(h3("Correlated Instrumental Variables"),
                   sliderInput(inputId = "rho",
                               label = "Correlation between dependend IVs (&rho;)",
                               min = 0, max = 1, value = 0.5, step = 0.1),
                   h3("Effects of Correlated IVs on X"),
                   numericInput(inputId = "depIV1_on_X",
                                label = HTML(paste0("Z<sub>", as.numeric(input$n_indIVs) + 1, "</sub> (&gamma;<sub>D1</sub>)")),
                                value = 0.1),
                   numericInput(inputId = "depIV2_on_X",
                                label = HTML(paste0("Z<sub>", as.numeric(input$n_indIVs) + 2, "</sub> (&gamma;<sub>D2</sub>)")),
                                value = 0.1)
      )

      return(outs)
    }

  })

  RVs <- reactiveValues()

  observeEvent(input$simulate_data, {
    updateTabsetPanel(
      inputId = "two_sample_bounds",
      selected = "Simulation Results"
    )

    if(input$n_indIVs > 0){
      RVs$indIVs_on_X <- sapply(1:input$n_indIVs,
                                function(i) input[[paste0("Z", i, "_on_X")]])
      RVs$indIVs_on_Y <- sapply(1:input$n_indIVs,
                                function(i) input[[paste0("Z", i, "_on_Y")]])
    } else {
      RVs$indIVs_on_X <- NULL
      RVs$indIVs_on_Y <- NULL
    }

    if(!input$invalid_IVs)
      RVs$indIVs_on_Y <- NULL

    # RVs$IVs_ps <- lapply(1:(as.numeric(input$n_cats)-1),
    #                      function(i) input[[paste0("IVs_ps", i)]])

    ## Create list with one element per ind. IV. The i'th element should be
    ## a vector of length input$n_cats with P(Z_i = j) as the j'th entry.
    RVs$IVs_ps <- map(1:as.numeric(input$n_indIVs),
                      function(x){
                        tmp <- lapply(0:(as.numeric(input$n_cats)-1),
                                      function(i){
                                        input[[paste0("IVs_ps", x, "_", i)]]
                                      }) %>% unlist()
                        return(c(tmp, 1-sum(tmp)))
                      })

    #print(paste("Before:", paste(RVs$IVs_ps, collapse = ", ")))

    # RVs$IVs_ps <- lapply(1:(as.numeric(input$n_cats)-1),
    #        function(i){
    #          lapply(1:as.numeric(input$n_indIVs),
    #                 function(i,k){
    #                   c(input[[paste0("IVs_ps", i, "_", k)]],
    #                     1 - sum(input[[paste0("IVs_ps", i, "_", k)]]))
    #                 })
    #        })

    #RVs$IVs_ps[length(RVs$IVs_ps) + 1] <- 1 - sum(RVs$IVs_ps)

    #print(paste("After:", paste(RVs$IVs_ps, collapse = ", ")))

    if(input$include_dep_IVs){
      RVs$depIVs_on_X <- c(input$depIV1_on_X, input$depIV2_on_X)
    } else {
      RVs$depIVs_on_X <- NULL
    }



    RVs$sim_data <- ATEBounds::simulate_data(
      sample_size = input$sample_size,
      indIVs_on_X = RVs$indIVs_on_X,
      indIVs_on_Y = RVs$indIVs_on_Y,
      IVs_ps = RVs$IVs_ps,
      U_on_X = input$U_on_X,
      U_on_Y = input$U_on_Y,
      X_on_Y = input$X_on_Y,
      rho = input$rho,
      depIVs_on_X = RVs$depIVs_on_X,
      dep_probs = ifelse(input$n_cats == 2, rep(0.5, 2), c(0.25, 0.5, 0.25)),
      X_intercept = -sum(RVs$indIVs_on_X),
      Y_intercept = -input$X_on_Y,
      U = distributions3::Normal()
    )

    readr::write_rds(RVs$sim_data, file = "sim_data_from_shiny.Rds")
  })

  observeEvent(RVs$sim_data, {

    RVs$results <- RVs$sim_data$simulated_data %>%
      pivot_longer(cols = starts_with("Z"),
                   names_to = "p", values_to = "Z") %>%
      nest_by(p) %>%
      summarize(probs = list(map(probs_from_data(dat = data,
                                                   X = X, Y = Y, Z = Z, data_format = "bivariate"),
                                 round, digits = 5)),
                .groups = "drop") %>%
      unnest_wider(probs) %>%
      mutate(p = str_remove(p, "Z") %>% as.numeric(),
             bounds = map2(gammas, thetas,
                           get_bounds, stop = FALSE, warning = FALSE),
             constraints_violated = map_lgl(bounds, "constraints_violated"),
             ints = map(bounds, "interval")) %>%
      select(-bounds) %>%
      unnest_wider(ints)


    output$tidy_results <- DT::renderDT({
      RVs$results %>%
        rowwise() %>%
        mutate(gammas = paste(gammas, collapse = ", "),
               thetas = paste(thetas, collapse = ", ")) %>%
        rename("P(X = 1 | Z = z)" = thetas,
               "P(Y = 1 | Z = z)" = gammas,
               "Constraints Violated" = constraints_violated) %>%
        arrange(p) %>%
        mutate(p = paste0("Z", p)) %>%
        rename(IV = p) %>%
        DT::datatable(
          # filter = "top",
          options = list(dom = "tip"),
          caption =
              paste("Note: P(X = 1 | Z = z) and P(Y = 1 | Z = z) gives values for z =",
                    paste(seq_along(RVs$results$gammas[[1]])-1, collapse = ", "),
                    "in that order.")
        )
    })

    output$simulated_data <- DT::renderDT({
      DT::datatable(RVs$sim_data$simulated_data,
                    filter = "top")
    })

    output$DAG_plotted <- renderPlot({
      # tidyDAG <- from_simulated_to_dag(RVs$sim_data)
      #
      # ggdag(tidyDAG) +
      #   geom_label_repel(aes(label = coef,
      #                        x = label_x,
      #                        y = label_y)) +
      #   theme_dag_blank() +
      #   ggtitle("DAG structure used for simulation")
      plot_DAG(RVs$sim_data)
    })

    output$plot_results <- renderPlot({
      plot_bounds(RVs$results, RVs$sim_data$coefficients)
    })

  })


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
    }
  })
}
