library(tidyverse)
library(flextable)
library(tinytable)
library(table1)

table_e1_synthetic = read_csv("../python/synthetic_null_generation/figures/table_e1-synthetic_data.csv")

flextable(table_e1_synthetic) |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-synthetic_data.png")

table_e1_oracle = table_e1_synthetic |>
  filter(str_detect(Dataset, "original"))
flextable(table_e1_oracle) |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-synthetic_data-oracle.png")

table_e1_ss = table_e1_synthetic |>
  filter(str_detect(Dataset, "Strimmer")) |>
  mutate(corr_type = "Schaefer-Strimmer", .after = Dataset)
flextable(table_e1_ss) |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-synthetic_data-schaefer_strimmer.png")

table_e1_empirical = table_e1_synthetic |>
  filter(str_detect(Dataset, "empirical")) |>
  mutate(corr_type = "Empirical", .after = Dataset)
flextable(table_e1_empirical) |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-synthetic_data-empirical.png")

# https://stackoverflow.com/questions/71035743/r-is-there-a-way-to-round-to-1-decimal-place-in-table1-for-continuous-variables
my.render.cont <- function(x) {
  with(stats.default(x), 
       c("",
         
         "Mean (SD)" = sprintf("%s (%s)",
                               round_pad(MEAN, 1),
                               round_pad(SD, 1)),
         
         "Median (Min, Max)" = sprintf("%s (%s, %s)",
                                       round_pad(MEDIAN, 1), 
                                       round_pad(MIN, 1), 
                                       round_pad(MAX, 1)))
  )
}

table_e1_only_synthetic = table_e1_ss |>
  bind_rows(table_e1_empirical)
t1_means_vars = table1(~ Genes + `Minimum seq. depth` + `Maximum seq. depth` + 
         `Median seq. depth` + `Zero counts (percentage)` + `Maximum count` + 
         `95% quantile` + `99% quantile` | corr_type, table_e1_only_synthetic,
       render.continuous = my.render.cont) 
t1_means_vars |>
  t1flex() |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-synthetic_data-mean_var.png")

table_e1_only_synthetic_means = table_e1_only_synthetic |>
  select(-c(Dataset)) |>
  group_by(corr_type) |>
  summarise(across(everything(), mean)) |>
  mutate(Dataset = paste("Mean of", corr_type, "estimate", sep = " "))

table_e1_oracle_with_means = bind_rows(table_e1_oracle, table_e1_only_synthetic_means)
table_e1_oracle_with_means |>
  flextable() |>
  save_as_image("../../seminar_paper-custom/figures/table_e1-oracle_with_means.png")
  