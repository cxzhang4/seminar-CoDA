library(tidyverse)
library(flextable)

table_e1_synthetic = read_csv("../python/reproduce_results/table_e1.csv")
flextable(table_e1_synthetic) |>
    save_as_image("../../seminar_paper-custom/figures/table_e1-reproduced.png")