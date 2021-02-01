latex_lower_bounds <- ACEBounds:::matrices_from_polymake %>%
  filter(n_z_levels == 3, !x_monotone, !y_monotone, data_format == "bivariate") %>%
  pull(matrix) %>% .[[1]] %>%
  filter(alpha == 1) %>%
  mutate(across(everything(), ~-.x),
         id = row_number()) %>%
  pivot_longer(cols = -c(alpha, id), values_to = "coef") %>%
  mutate(Z = paste("Z =", str_sub(name, start = -1)),
         XY_value = str_sub(name, start = -2, end = -2),
         XY = if_else(str_detect(name, "gamma"), "Y = 1", "X = 1")) %>%
  select(-name, -alpha) %>%
  pivot_wider(names_from = XY_value, values_from = coef) %>%
  mutate(coef = `1` - `0`,
         p = case_when(coef == 1 ~ paste0("P(", XY, " | ", Z, ")"),
                       coef == -1 ~ paste0("-P(", XY, " | ", Z, ")"),
                       coef == 0 ~ "",
                       TRUE ~ paste0(coef, "\\cdot P(", XY, " | ", Z, ")"))) %>%
  rename(const = `0`) %>%
  group_by(id) %>%
  mutate(p = if_else(row_number() != 1 & str_sub(p, end = 1) != "-" & p != "",
                     paste(" + ", p), p)) %>%
  summarize(p = paste(p, collapse = ""),
            const = sum(const)) %>%
  mutate(out = case_when(const < 0 ~ paste(p, const),
                         const > 0 ~ paste(p, const, sep = " + "),
                         TRUE ~ p),
         out = if_else(str_sub(out, start = 1, end = 3) == " + ",
                       str_sub(out, start = 4),
                       out)) %>%
  pull(out)

latex_upper_bounds <- ACEBounds:::matrices_from_polymake %>%
  filter(n_z_levels == 3, !x_monotone, !y_monotone, data_format == "bivariate") %>%
  pull(matrix) %>% .[[1]] %>%
  filter(alpha == -1) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = -c(alpha, id), values_to = "coef") %>%
  mutate(Z = paste("Z =", str_sub(name, start = -1)),
         XY_value = str_sub(name, start = -2, end = -2),
         XY = if_else(str_detect(name, "gamma"), "Y = 1", "X = 1")) %>%
  select(-name, -alpha) %>%
  pivot_wider(names_from = XY_value, values_from = coef) %>%
  mutate(coef = `1` - `0`,
         p = case_when(coef == 1 ~ paste0("P(", XY, " | ", Z, ")"),
                       coef == -1 ~ paste0("-P(", XY, " | ", Z, ")"),
                       coef == 0 ~ "",
                       TRUE ~ paste0(coef, "\\cdot P(", XY, " | ", Z, ")"))) %>%
  rename(const = `0`) %>%
  group_by(id) %>%
  mutate(p = if_else(row_number() != 1 & str_sub(p, end = 1) != "-" & p != "",
                     paste(" + ", p), p)) %>%
  summarize(p = paste(p, collapse = ""),
            const = sum(const)) %>%
  mutate(out = case_when(const < 0 ~ paste(p, const),
                         const > 0 ~ paste(p, const, sep = " + "),
                         TRUE ~ p),
         out = if_else(str_sub(out, start = 1, end = 3) == " + ",
                       str_sub(out, start = 4),
                       out)) %>%
  pull(out) %>%
  str_trim()


latex_constraints <- ACEBounds:::matrices_from_polymake %>%
  filter(n_z_levels == 3, !x_monotone, !y_monotone, data_format == "bivariate") %>%
  pull(matrix) %>% .[[1]] %>%
  filter(alpha == 0) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = -c(alpha, id), values_to = "coef") %>%
  mutate(Z = paste("Z =", str_sub(name, start = -1)),
         XY_value = str_sub(name, start = -2, end = -2),
         XY = if_else(str_detect(name, "gamma"), "Y = 1", "X = 1")) %>%
  select(-name, -alpha) %>%
  pivot_wider(names_from = XY_value, values_from = coef) %>%
  mutate(coef = `1` - `0`,
         p = case_when(coef == 1 ~ paste0("P(", XY, " | ", Z, ")"),
                       coef == -1 ~ paste0("-P(", XY, " | ", Z, ")"),
                       coef == 0 ~ "",
                       TRUE ~ paste0(coef, "\\cdot P(", XY, " | ", Z, ")"))) %>%
  rename(const = `0`) %>%
  group_by(id) %>%
  mutate(p = if_else(row_number() != 1 & str_sub(p, end = 1) != "-" & p != "",
                     paste(" + ", p), p)) %>%
  summarize(p = paste(p, collapse = ""),
            const = sum(const)) %>%
  mutate(out = case_when(const < 0 ~ paste(p, const),
                         const > 0 ~ paste(p, const, sep = " + "),
                         TRUE ~ p),
         out = if_else(str_sub(out, start = 1, end = 3) == " + ",
                       str_sub(out, start = 4),
                       out)) %>%
  pull(out) %>% str_trim %>% sort()
