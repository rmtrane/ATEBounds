library(tidyverse)

epi_pdf <- here::here("epidemiology/epi.pdf")

lines <- pdftools::pdf_text(epi_pdf) %>% read_lines() %>% str_trim()

references_start <- str_detect(lines, "^References")
sdc_start <- str_detect(lines, "^A [[:space:]]+ Supplemental")

lines_w_pagenumbers <- tibble(lines = lines) %>%
  mutate(page_number = as.numeric(str_extract(lines, pattern = "^[1-9][0-9]*$"))) %>%
  fill(page_number, .direction = "up")

main_ends <- lines_w_pagenumbers %>%
  filter(lines != "") %>%
  filter(row_number() < which(str_detect(lines, pattern = "^A [[:space:]]+ Supplemental")),
         page_number > 1) %>%
  tail(n=1) %>% pull(page_number)

pdftools::pdf_subset(input = epi_pdf, pages = 1:main_ends,
                     output = here::here("epidemiology/epi_main.pdf"))

pdftools::pdf_subset(input = epi_pdf, pages = (main_ends + 1):pdftools::pdf_length(epi_pdf),
                     output = here::here("epidemiology/epi_sdc.pdf"))


### SIM
sim_pdf <- here::here("sim_submission/sim_submission.pdf")

lines <- pdftools::pdf_text(sim_pdf) %>% read_lines() %>% str_trim()

references_start <- str_detect(lines, "^References")
appendix_start <- str_detect(lines, "^Appendix$")

lines_w_pagenumbers <- tibble(lines = lines) %>%
  mutate(page_number_line = str_detect(lines, pattern = "^[0-9]+ [[:space:]]+ Trane and Kang$|^Trane and Kang [[:space:]]+ [0-9]+$|^[0-9]+$"),
         page_number = if_else(page_number_line, str_remove(lines, "Trane and Kang") %>% str_trim() %>% as.numeric(),
                               NA_real_)) %>%
  fill(page_number, .direction = "down") %>%
  mutate(page_number = if_else(is.na(page_number), 1, page_number) %>% cummax())

main_ends <- lines_w_pagenumbers %>%
  filter(lines != "") %>%
  filter(row_number() < which(str_detect(tolower(lines), pattern = "^appendix$")),
         page_number > 1) %>%
  tail(n=1) %>% pull(page_number)

pdftools::pdf_subset(input = sim_pdf, pages = 1:main_ends,
                     output = here::here("sim_submission/sim_main.pdf"))

pdftools::pdf_subset(input = sim_pdf, pages = (main_ends + 1):pdftools::pdf_length(sim_pdf),
                     output = here::here("sim_submission/sim_appendix.pdf"))
