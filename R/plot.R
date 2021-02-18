#' Create tidy DAG from simulation
#'
from_simulated_to_dag <- function(simulated_data, layout = "nicely"){
  if(!"ggdag" %in% rownames(installed.packages()))
    stop("You have to install the ggdag package to use this function.")

  coefs <- simulated_data$coefficients

  tmp <- coefs %>%
    filter(str_detect(effect, "_on_")) %>%
    mutate(edges = str_replace_all(effect, pattern = "_on_", "->") %>% paste0(., ";"))

  all_nodes <- tmp$effect %>% str_split("_on_", simplify = T) %>% .[,1] %>% unique()

  DAG <- paste0("dag{",
                paste(tmp$edges, collapse = ""),
                "}") %>%
    dagitty::dagitty()

  tidy_dag <- ggdag::tidy_dagitty(DAG,
                                  layout = layout)

  # Check if any IVs are placed at different level than the rest
  # if(n_distinct(filter(tidy_dag$data, str_starts(name, "Z"))$y) > 1){
  #   # y's
  #   ys <- filter(tidy_dag$data, str_starts(name, "Z"))$y
  #   # unique y's
  #   uniq_ys <- unique(ys)
  #   # most common y
  #   Z_y <- uniq_ys[which.max(tabulate(match(ys, uniq_ys)))]
  #
  #   tidy_dag$data <- tidy_dag$data %>%
  #     mutate(y = if_else(str_starts(name, "Z"), Z_y, y))
  # }

  tidy_dag$data <- tidy_dag$data %>%
    left_join(
      coefs %>%
        separate(effect, into = c("name", "to"), sep = "_on_", fill = "right") %>%
        filter(!is.na(to))
    )

  any_corrZs <- str_detect(coefs$effect,
                           pattern = "Z[0-9]_on_Z[0-9]")

  if(any(any_corrZs)){
    corrZs <- str_split(coefs$effect[any_corrZs], pattern = "_on_")[[1]]

    tidy_dag$data <- tidy_dag$data %>%
      mutate(direction = if_else(name %in% corrZs & to %in% corrZs,
                                 factor("<->", levels = levels(direction)),
                                 direction))
  }

  tidy_dag$data <- tidy_dag$data %>%
    mutate(label_x = (x + xend)/2,
           label_y = (y + yend)/2)

  return(tidy_dag)
}

#' Plot DAG from simulation
#'
#' @export
plot_DAG <- function(sim_data){

  if(!"ggdag" %in% rownames(installed.packages()))
    stop("You have to install the ggdag package to use this function.")

  tidyDAG <- from_simulated_to_dag(sim_data)

  ggdag::ggdag(tidyDAG) +
    # ggdag::geom_dag_label_repel(aes(label = coef,
    #                                 x = label_x,
    #                                 y = label_y),
    #                             force = 0) +
    geom_label(data = filter(tidyDAG, !is.na(coef)),
               aes(label = coef,
                   x = label_x,
                   y = label_y)) +
    ggdag::theme_dag_blank() +
    ggtitle("DAG structure used for simulation")
}
