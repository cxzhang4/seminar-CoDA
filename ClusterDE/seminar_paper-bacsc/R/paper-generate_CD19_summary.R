library(tidyverse)
library(flextable)

CD19 = read_csv("../../seminar_paper-custom/CD19.csv")
CD19 |>
  count(Species) |>
  flextable() |>
  save_as_image("../../seminar_paper-custom/figures/CD19_species_counts.png")
  
CD19 |>
  head(10) |>
  flextable() |>
  save_as_image("../../seminar_paper-custom/figures/CD19_head.png")
