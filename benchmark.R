library(bench)
library(tidyverse)

# -------------------------------------------------------
# Benchmark: vary number of items (p)
# -------------------------------------------------------

results_p <- bind_rows(
  lapply(c(20, 50, 100, 200, 500), function(p) {
    lapply(c("kw", "permute"), function(method) {
      gc()
      
      result <- bench::mark(
        {
          data <- treatment_data_generation_3d(p, 20, 0.8, 0, 0, 1:10)
          VSAT(data, method = method, n_perm = 200)
        },
        iterations = 5,
        check      = FALSE,
        filter_gc  = FALSE
      )
      gc()
      result$p      <- p
      result$method <- method
      result
    }) |> bind_rows()
  }) |> bind_rows()
)

saveRDS(results_p, "eval_results/item_time_benchmark.RDS")

# -------------------------------------------------------
# Benchmark: vary number of samples (n)
# -------------------------------------------------------

results_n <- bind_rows(
  lapply(c(10, 20, 36, 50, 100), function(n) {
    lapply(c("kw", "permute"), function(method) {
      gc()
      
      result <- bench::mark(
        {
          data <- treatment_data_generation_3d(100, n, 0.8, 0, 0, 1:10)
          VSAT(data, method = method, n_perm = 200)
        },
        iterations = 5,
        check      = FALSE,
        filter_gc  = FALSE
      )
      gc()
      result$n      <- n
      result$method <- method
      result
    }) |> bind_rows()
  }) |> bind_rows()
)

saveRDS(results_n, "eval_results/sample_time_benchmark.RDS")

# -------------------------------------------------------
# Benchmark: vary number of permutations (permute only)
# -------------------------------------------------------

results_perm <- bind_rows(
  lapply(c(100, 200, 500, 1000), function(n_perm) {
    gc()
    iters <- if (n_perm >= 500) 3L else 5L
    
    result <- bench::mark(
      {
        data <- treatment_data_generation_3d(100, 20, 0.8, 0, 0, 1:10)
        VSAT(data, method = "permute", n_perm = n_perm)
      },
      iterations = 5,
      check      = FALSE,
      filter_gc  = FALSE
    )
    gc()
    result$n_perm <- n_perm
    result
  }) |> bind_rows()
)

saveRDS(results_perm, "eval_results/num_perm_benchmark.RDS")

# -------------------------------------------------------
# Benchmark: informed starting point (3/10 true DC taxa)
# -------------------------------------------------------

results_A0 <- bind_rows(
  lapply(c(20, 50, 100, 200, 500), function(p) {
    lapply(c("kw", "permute"), function(method) {
      gc()
      
      result <- bench::mark(
        {
          data <- treatment_data_generation_3d(p, 20, 0.8, 0, 0, 1:10)
          VSAT(data, A0 = c("V1", "V3", "V5"), method = method, n_perm = 200)
        },
        iterations = 5,
        check      = FALSE,
        filter_gc  = FALSE
      )
      gc()
      result$p      <- p
      result$method <- method
      result
    }) |> bind_rows()
  }) |> bind_rows()
)

saveRDS(results_A0, "eval_results/informed_start_benchmark.RDS")

# -------------------------------------------------------
# Load results (if re-running plots without re-benchmarking)
# -------------------------------------------------------

results_p    <- readRDS("eval_results/item_time_benchmark.RDS")
results_n    <- readRDS("eval_results/sample_time_benchmark.RDS")
results_perm <- readRDS("eval_results/num_perm_benchmark.RDS")
results_A0   <- readRDS("eval_results/informed_start_benchmark.RDS")

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
    subtitle = expression("n = 36, " ~ rho[1] == 0.8 ~ ", n_perm = 200"),
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
    subtitle = expression("p = 100, " ~ rho[1] == 0.8 ~ ", n_perm = 200"),
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
    subtitle = expression("p = 100, n = 20, " ~ rho[1] == 0.8),
    x        = "Number of Permutations",
    y        = "Median Runtime (seconds)"
  ) +
  vsat_theme()

# greedy vs informed A0 - both methods
bind_rows(
  results_p  |> mutate(init = "Greedy"),
  results_A0 |> mutate(init = "Informed Start")
) |>
  mutate(
    median_time = as.double(median),
    method      = case_when(
      method == "kw"      ~ "Kruskal-Wallis",
      method == "permute" ~ "Permutation"
    ),
    init = factor(init, levels = c("Greedy", "Informed Start"))
  ) |>
  ggplot(aes(x = p, y = median_time,
             color = method, linetype = init,
             group = interaction(method, init))) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.9) +
  vsat_colors +
  scale_linetype_manual(
    values = c("Greedy"                      = "solid",
               "Informed Start" = "dashed")
  ) +
  labs(
    title    = "Runtime: Greedy vs Informed Initialization",
    subtitle = expression("n = 20, " ~ rho[1] == 0.8 ~ ", n_perm = 200"),
    x        = "Number of Items (p)",
    y        = "Median Runtime (seconds)",
    color    = "Method",
    linetype = "Initialization"
  ) +
  vsat_theme() +
  theme(
    legend.position = "right"
  )

