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
