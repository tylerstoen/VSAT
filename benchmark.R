library(bench)
library(tidyverse)

# -------------------------------------------------------
# Benchmark: vary number of items (p)
# -------------------------------------------------------

results_p <- bench::press(
  p = c(20, 50, 100, 200, 500),
  method = c("kw", "permute"),
  {
    bench::mark(
      {
        data <- treatment_data_generation_3d(p, 20, 0.8, 0, 0, 1:10)
        VSAT(data, method = method, n_perm = 200)
      },
      iterations = 5,
      check      = FALSE
    )
  }
)

saveRDS(results_p, "eval_results/item_time_benchmark.RDS")

results_p <- readRDS("eval_results/item_time_benchmark.RDS")

# -------------------------------------------------------
# Benchmark: vary number of samples (n)
# -------------------------------------------------------

results_n <- bench::press(
  n = c(10, 20, 36, 50, 100),
  method = c("kw", "permute"),
  {
    bench::mark(
      {
        data <- treatment_data_generation_3d(100, n, 0.8, 0, 0, 1:10)
        VSAT(data, method = method, n_perm = 200)
      },
      iterations = 5,
      check      = FALSE
    )
  }
)

saveRDS(results_n, "eval_results/sample_time_benchmark.RDS")

results_n <- readRDS("eval_results/sample_time_benchmark.RDS")

# -------------------------------------------------------
# Benchmark: vary number of permutations (permute method only)
# -------------------------------------------------------

results_perm <- bench::press(
  n_perm = c(100, 200, 500, 1000),
  {
    bench::mark(
      {
        data <- treatment_data_generation_3d(100, 20, 0.8, 0, 0, 1:10)
        VSAT(data, method = "permute", n_perm = n_perm)
      },
      iterations = 5,
      check      = FALSE
    )
  }
)

saveRDS(results_perm, "eval_results/num_perm_benchmark.RDS")

# -------------------------------------------------------
# Plots
# -------------------------------------------------------

# items (p) - both methods
results_p |>
  mutate(
    median_time = as.double(median),
    method      = case_when(
      method == "kw"      ~ "Kruskal-Wallis",
      method == "permute" ~ "Permutation"
    )
  ) |>
  ggplot(aes(x = p, y = median_time, color = method, group = method)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  vsat_colors +
  labs(
    title    = "Runtime by Number of Items",
    subtitle = expression("n = 36, " ~ rho[1] ~ " = 0.8, n_perm = 200"),
    x        = "Number of Items (p)",
    y        = "Median Runtime (seconds)",
    color    = "Method"
  ) +
  vsat_theme()

# samples (n) - both methods
results_n |>
  mutate(
    median_time = as.double(median),
    method      = case_when(
      method == "kw"      ~ "Kruskal-Wallis",
      method == "permute" ~ "Permutation"
    )
  ) |>
  ggplot(aes(x = n, y = median_time, color = method, group = method)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  vsat_colors +
  labs(
    title    = "Runtime by Sample Size",
    subtitle = "p = 100," ~ rho[1] ~ "= 0.8, n_perm = 200",
    x        = "Samples per Group (n)",
    y        = "Median Runtime (seconds)",
    color    = "Method"
  ) +
  vsat_theme()

# permutations - permute method only
results_perm |>
  mutate(median_time = as.double(median)) |>
  ggplot(aes(x = n_perm, y = median_time)) +
  geom_point(size = 2.5, color = "#619CFF") +
  geom_line(linewidth = 0.9, color = "#619CFF") +
  labs(
    title    = "Permutation Method Runtime by Number of Permutations",
    subtitle = "p = 100, n = 36, ρ₁ = 0.8, ρ₂ = ρ₃ = 0",
    x        = "Number of Permutations",
    y        = "Median Runtime (seconds)"
  ) +
  vsat_theme()