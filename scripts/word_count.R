library(pdftools)
library(tidyverse)

pdf_word_count <- function(pdf_file = here::here("epidemiology/epi.pdf"),
                           main_starts = "^1 [[:space:]]+ Introduction",
                           appendix_starts = "^A [[:space:]]+ Appendix"){

  doc_lines <- pdf_text(pdf_file) %>%
    readr::read_lines() %>%
    str_trim()

  no_appendix <- tibble(lines = doc_lines) %>%
    filter(str_detect(lines, "[A-Za-z]"),
           !str_detect(lines, "/Users/"))

  oneliner <- no_appendix %>%
    pull(lines) %>%
    paste(collapse = " ")

  str_count(oneliner, "\\w+")
}

pdf_word_count()

words_per_line <- tibble(lines = doc_lines) %>%
  filter(str_detect(lines, "[A-Za-z]")) %>%
  mutate(section_starts = str_detect(lines, "^([0-9\\.]+|A[[[:space:]]\\.0-9]+) [[:space:]] [[:space:]]+"), # (Appendix|Introduction|Methods|Review|Notation|Properties|Using|Method|Discussion|Results|IV|Length|Would|Characterizing)"),
         section = if_else(section_starts, str_trim(str_extract(lines, "[0-9\\.]+|^A[0-9\\. [[:space:]]]+")), NA_character_)) %>%
  fill(section) %>%
  mutate(section = if_else(is.na(section), "0", section)) %>%
  mutate(n_words = str_count(lines, "\\w+"))

words_per_line %>% View()

words_per_line %>%
  group_by(section) %>%
  summarize(n_words = sum(n_words)) %>%
  print(n = Inf)


# Without abstract, section headers, and large equations:
words_per_line %>%
  filter(!row_number() %in% c(110:120)) %>%
  group_by(section) %>%
  summarize(n_words = sum(n_words)) %>%
  print(n = Inf) %>% ungroup() %>%
  filter(row_number() > which(str_trim(section) == "0"),
         row_number() < which(str_trim(section) == "A")) %>%
  summarize(total = sum(n_words))

# oneliner <- epi_pdf_text %>%
#   pull(lines) %>%
#   paste(collapse = " ")
#
# str_count(oneliner, "[\\w\']+")
#
# str_count(tmp, "\\w+")
#
# fifty <- "A wonderful serenity has taken possession of my entire soul, like these sweet mornings of spring which I enjoy with my whole heart. I am alone, and feel the charm of existence in this spot, which was created for the bliss of souls like mine. I am so happy, my"
# fivehundred <- "A wonderful serenity has taken possession of my entire soul, like these sweet mornings of spring which I enjoy with my whole heart. I am alone, and feel the charm of existence in this spot, which was created for the bliss of souls like mine. I am so happy, my dear friend, so absorbed in the exquisite sense of mere tranquil existence, that I neglect my talents. I should be incapable of drawing a single stroke at the present moment; and yet I feel that I never was a greater artist than now. When, while the lovely valley teems with vapour around me, and the meridian sun strikes the upper surface of the impenetrable foliage of my trees, and but a few stray gleams steal into the inner sanctuary, I throw myself down among the tall grass by the trickling stream; and, as I lie close to the earth, a thousand unknown plants are noticed by me: when I hear the buzz of the little world among the stalks, and grow familiar with the countless indescribable forms of the insects and flies, then I feel the presence of the Almighty, who formed us in his own image, and the breath of that universal love which bears and sustains us, as it floats around us in an eternity of bliss; and then, my friend, when darkness overspreads my eyes, and heaven and earth seem to dwell in my soul and absorb its power, like the form of a beloved mistress, then I often think with longing, Oh, would I could describe these conceptions, could impress upon paper all that is living so full and warm within me, that it might be the mirror of my soul, as my soul is the mirror of the infinite God! O my friend -- but it is too much for my strength -- I sink under the weight of the splendour of these visions! A wonderful serenity has taken possession of my entire soul, like these sweet mornings of spring which I enjoy with my whole heart. I am alone, and feel the charm of existence in this spot, which was created for the bliss of souls like mine. I am so happy, my dear friend, so absorbed in the exquisite sense of mere tranquil existence, that I neglect my talents. I should be incapable of drawing a single stroke at the present moment; and yet I feel that I never was a greater artist than now. When, while the lovely valley teems with vapour around me, and the meridian sun strikes the upper surface of the impenetrable foliage of my trees, and but a few stray gleams steal into the inner sanctuary, I throw myself down among the tall grass by the trickling stream; and, as I lie close to the earth, a thousand unknown plants are noticed by me: when I hear the buzz of the little world among the stalks, and grow familiar with the"
#
# str_count(fifty, " ")
# str_count(fifty, "\\w+")
#
# str_count(fivehundred, "\\w+")
# str_count(fivehundred, "\\w+")