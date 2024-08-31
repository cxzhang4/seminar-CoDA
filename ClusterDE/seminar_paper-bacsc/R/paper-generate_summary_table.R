library(tidyverse)
library(flextable)
library(tinytable)

table_e1_synthetic = read_csv("../python/synthetic_null_generation/figures/table_e1-synthetic_data.csv")

flextable(table_e1_synthetic) |>
  save_as_image("figures/table_e1-synthetic_data.png")

table_e1_ss = table_e1_synthetic |>
  filter(str_detect(Dataset, "Strimmer"))

table_e1_empirical = table_e1_synthetic |>
  filter(str_detect(Dataset, "empirical"))
