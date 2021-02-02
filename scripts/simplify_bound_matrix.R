library(tidyverse)

simplify_matrix <- function(poly_matrix, upper_or_lower = c("upper", "lower")){
  poly_matrix %>%
    select(where(~sum(abs(.x)) > 0)) %>%
    filter(alpha == (1 - 2*as.numeric(upper_or_lower == "lower"))) %>%
    select(-alpha) %>%
    mutate(across(everything(), ~(1-2*as.numeric(upper_or_lower == "lower"))*.x),
           id = row_number()) %>%
    pivot_longer(cols = -id) %>%
    mutate(z = as.numeric(str_sub(name, start = -1, end = -1)) - 1,
           xy = as.numeric(str_sub(name, start = -2, end = -2)),
           name = if_else(str_sub(name, start = 1, end = -3) == "gamma", "Y", "X")) %>%
    pivot_wider(names_from = xy, values_from = value, values_fill = 0) %>%
    mutate(`1` = `1` - `0`,
           const = `0`) %>%
    select(-`0`) %>%
    filter(abs(`1`) + abs(const) > 0) %>%
    group_by(id) %>%
    mutate(`1` = case_when(`1` == -1 & row_number() == 1 ~ "-",
                           `1` == -1 & row_number() > 1 ~ "- ",
                           `1` == 1 & row_number() == 1 ~ "",
                           `1` == 1 & row_number() > 1 ~ "+ ",
                           row_number() > 1 & sign(`1`) == 1 ~ paste("+", `1`),
                           `1` == 0 ~ "",
                           TRUE ~ as.character(`1`)),
           terms = glue::glue("{`1`}P({name} = 1 | Z = {z})")) %>%
    summarize(const = case_when(sum(const) > 0 ~ paste("+", sum(const)),
                                sum(const) < 0 ~ paste("-", abs(sum(const))),
                                TRUE ~ ""),
              expression = paste(paste(terms, collapse = " "), const)) %>%
    select(id, expression)
}

ATEBounds:::matrices_from_polymake %>%
  filter(x_monotone, y_monotone, data_format == "bivariate") %>%
  mutate(simplified_matrix_lower = map(matrix, simplify_matrix, upper_or_lower = "lower"),
         simplified_matrix_upper = map(matrix, simplify_matrix, upper_or_lower = "upper")) %>%
  pull(simplified_matrix_upper)
