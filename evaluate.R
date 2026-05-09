library(MASS)
library(tidyverse)
library(furrr)


# Scenario 1: single signal, varying strength
scenario_1 <- crossing(
  scenario = "Single Signal",
  p        = 100,
  n        = 36,
  set_size = 10,
  rho_1    = c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9),
  rho_2    = 0,
  rho_3    = 0,
  method   = c("kw", "permute"),
  rep      = 1:50
)

# Scenario 2: gradient signal
scenario_2 <- crossing(
  scenario = "Signal Gradient",
  p        = 100,
  n        = 36,
  set_size = 10,
  rho_1    = c(0.9, 0.8, 0.7, 0.6, 0.5),   # strongest signal
  rho_2    = NA,                              # placeholder
  rho_3    = NA,
  method   = c("kw", "permute"),
  rep      = 1:50
) |>
  mutate(
    rho_2 = rho_1 * 0.5,    # half the signal of rho_1
    rho_3 = rho_1 * 0.25    # quarter the signal of rho_1
  )

scenario_3 <- crossing(
  scenario = "Signal Gradient Constant",
  p        = 100,
  n        = 36,
  set_size = 10,
  rho_1    = c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9),  # varies
  rho_2    = 0.3,                                 # fixed background
  rho_3    = 0.1,                                 # fixed background
  method   = c("kw", "permute"),
  rep      = 1:50
)

# Scenario 4: varying p (population size)
scenario_4 <- crossing(
  scenario = "Population Size",
  p        = c(20, 50, 100, 250),
  n        = 36,
  set_size = 10,
  rho_1    = 0.6,
  rho_2    = 0,
  rho_3    = 0,
  method   = c("kw", "permute"),
  rep      = 1:50
)

# Scenario 5: varying n (sample size)
scenario_5 <- crossing(
  scenario = "Sample Size",
  p        = 100,
  n        = c(10, 20, 35, 50),
  set_size = 10,
  rho_1    = 0.6,
  rho_2    = 0,
  rho_3    = 0,
  method   = c("kw", "permute"),
  rep      = 1:50
)

# Scenario 6: varying set size
scenario_6 <- crossing(
  scenario = "Set Size",
  p        = 100,
  n        = 36,
  set_size = c(3, 5, 10, 15, 20),
  rho_1    = 0.6,
  rho_2    = 0,
  rho_3    = 0,
  method   = c("kw", "permute"),
  rep      = 1:50,
  
)

# Scenario null: no signal
scenario_null <- crossing(
  scenario = "Null",
  p        = 100,
  n        = 36,
  set_size = 10,
  rho_1    = seq(0, 0.9, by = 0.1),  # same as DCM x-axis
  rho_2    = seq(0, 0.9, by = 0.1),  # all equal - non-differential
  rho_3    = seq(0, 0.9, by = 0.1),
  method   = c("kw", "permute"),
  rep      = 1:25
) |>
  filter(rho_1 == rho_2 & rho_2 == rho_3)  # keep only equal rho rows




evaluate_vsat <- function(p, n, set_size, rho_1, rho_2, rho_3, 
                          method, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  true_set <- paste0("V", 1:set_size)
  
  data <- treatment_data_generation_3d(
    rows       = p,
    cols       = n,
    A_cov      = rho_1,   # explicit mapping here
    B_cov      = rho_2,
    C_cov      = rho_3,
    test_group = 1:set_size
  )
  
  result <- VSAT(data, method = method, n_perm = 200)
  
  if (!is.character(result) || length(result) == 0) {
    return(tibble(precision = NA_real_, recall = NA_real_,
                  f1 = NA_real_, empty = TRUE, returned_size = 0))
  }
  
  tp        <- length(intersect(result, true_set))
  precision <- tp / length(result)
  recall    <- tp / set_size
  f1        <- if (precision + recall == 0) 0 else
    2 * precision * recall / (precision + recall)
  
  tibble(
    precision     = precision,
    recall        = recall,
    f1            = f1,
    empty         = FALSE,
    returned_size = as.numeric(length(result))
  )
}

run_scenario <- function(scenario_grid, filename) {
  
  message("Starting: ", filename)
  
  results <- future_pmap(scenario_grid, function(...) {
    args <- list(...)
    
    # drop columns that aren't function arguments
    vsat_args <- args[c("p", "n", "set_size", "rho_1", "rho_2", "rho_3", "method")]
    vsat_args$seed <- args$rep
    
    res <- do.call(evaluate_vsat, vsat_args)
    bind_cols(as_tibble(args), res)
  }, .options = furrr_options(seed = TRUE)) |> list_rbind()
  
  saveRDS(results, filename)
  message("Saved: ", filename)
  return(results)
}

#results_s1 <- run_scenario(scenario_1, "results_single_signal.rds")
#results_s2 <- run_scenario(scenario_2, "results_gradient_signal.rds")
#results_s3 <- run_scenario(scenario_3, "results_gradient.rds")
#results_s4 <- run_scenario(scenario_4, "results_pop_size.rds")
#results_s5 <- run_scenario(scenario_5, "results_sample_size.rds")
#results_s6 <- run_scenario(scenario_6, "results_set_size.rds")
results_null <- run_scenario(scenario_null, "results_null.rds")


# -------------------------------------------------------
# Shared helper function
# -------------------------------------------------------

make_summary <- function(results, group_var) {
  results |>
    group_by({{ group_var }}, method) |>
    summarise(
      mean_recall    = mean(recall, na.rm = TRUE),
      mean_precision = mean(precision, na.rm = TRUE),
      mean_f1        = mean(f1, na.rm = TRUE),
      empty_rate     = mean(empty),
      .groups        = "drop"
    ) |>
    mutate(method = case_when(
      method == "kw"      ~ "Kruskal-Wallis",
      method == "permute" ~ "Permutation"
    )) |>
    pivot_longer(
      cols      = c(mean_recall, mean_precision, mean_f1, empty_rate),
      names_to  = "metric",
      values_to = "value"
    ) |>
    mutate(
      metric = case_when(
        metric == "mean_recall"    ~ "Recall",
        metric == "mean_precision" ~ "Precision",
        metric == "mean_f1"        ~ "F1",
        metric == "empty_rate"     ~ "Empty Set Rate"
      ),
      metric = factor(metric, levels = c(
        "Recall", "Precision", "F1", "Empty Set Rate"
      ))
    )
}

# shared plot theme
vsat_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# shared color scale
vsat_colors <- scale_color_manual(values = c(
  "Kruskal-Wallis" = "#F8766D",
  "Permutation"    = "#619CFF"
))

# shared FDR reference line — only appears in FDR panel
fdr_line <- geom_hline(
  data      = data.frame(metric = factor("False Discovery Rate",
                                         levels = c("True Positive Rate",
                                                    "False Discovery Rate",
                                                    "Empty Set Rate")),
                         yintercept = 0.05),
  aes(yintercept = yintercept),
  linetype  = "dashed",
  color     = "black",
  linewidth = 0.7
)

# -------------------------------------------------------
# Shared helper function
# -------------------------------------------------------

make_summary <- function(results, group_var) {
  results |>
    group_by({{ group_var }}, method) |>
    summarise(
      mean_recall    = mean(recall, na.rm = TRUE),
      mean_precision = mean(precision, na.rm = TRUE),
      mean_f1        = mean(f1, na.rm = TRUE),
      empty_rate     = mean(empty),
      .groups        = "drop"
    ) |>
    mutate(method = case_when(
      method == "kw"      ~ "Kruskal-Wallis",
      method == "permute" ~ "Permutation"
    )) |>
    pivot_longer(
      cols      = c(mean_recall, mean_precision, mean_f1, empty_rate),
      names_to  = "metric",
      values_to = "value"
    ) |>
    mutate(
      metric = case_when(
        metric == "mean_recall"    ~ "Recall",
        metric == "mean_precision" ~ "Precision",
        metric == "mean_f1"        ~ "F1",
        metric == "empty_rate"     ~ "Empty Set Rate"
      ),
      metric = factor(metric, levels = c(
        "Recall", "Precision", "F1", "Empty Set Rate"
      ))
    )
}

# shared plot theme
vsat_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# shared color scale
vsat_colors <- scale_color_manual(values = c(
  "Kruskal-Wallis" = "#F8766D",
  "Permutation"    = "#619CFF"
))

# -------------------------------------------------------
# Single Signal
# -------------------------------------------------------

summary_long_s1 <- make_summary(readRDS("results_single_signal.rds"), rho_1)

ggplot(summary_long_s1, aes(x = rho_1, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, ncol = 2) +
  vsat_colors +
  scale_x_continuous(breaks = unique(summary_long_s1$rho_1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "VSAT Performance Under Single Signal",
    subtitle = expression("p = 100, n = 36, set size = 10," ~
                            rho[2] == 0 ~ "," ~ rho[3] == 0),
    x        = expression("Signal strength" ~ rho[1]),
    y        = NULL,
    color    = "Method"
  ) +
  vsat_theme()

# -------------------------------------------------------
# Signal Gradient (proportional)
# -------------------------------------------------------

summary_long_s2 <- make_summary(readRDS("results_gradient_signal.rds"), rho_1)

ggplot(summary_long_s2, aes(x = rho_1, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, ncol = 2) +
  vsat_colors +
  scale_x_continuous(
    breaks = unique(summary_long_s2$rho_1),
    labels = function(x) parse(text = paste0("rho[1]==", x))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "VSAT Performance Under Proportional Signal Gradient",
    subtitle = expression("p = 100, n = 36, set size = 10," ~
                            rho[2] == frac(rho[1], 2) ~ "," ~
                            rho[3] == frac(rho[1], 4)),
    x        = expression("Primary signal strength" ~ rho[1]),
    y        = NULL,
    color    = "Method"
  ) +
  vsat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------------------------------------
# Signal Gradient (constant background)
# -------------------------------------------------------

summary_long_s3 <- make_summary(readRDS("results_gradient.rds"), rho_1)

ggplot(summary_long_s3, aes(x = rho_1, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, ncol = 2) +
  vsat_colors +
  scale_x_continuous(
    breaks = unique(summary_long_s3$rho_1),
    labels = function(x) parse(text = paste0("rho[1]==", x))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "VSAT Performance Under Constant Background Signal Gradient",
    subtitle = expression("p = 100, n = 36, set size = 10," ~
                            rho[2] == 0.3 ~ "," ~ rho[3] == 0.1),
    x        = expression("Primary signal strength" ~ rho[1]),
    y        = NULL,
    color    = "Method"
  ) +
  vsat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------------------------------------
# Item Population Size
# -------------------------------------------------------

summary_long_s4 <- make_summary(readRDS("results_pop_size.rds"), p)

ggplot(summary_long_s4, aes(x = p, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, ncol = 2) +
  vsat_colors +
  scale_x_continuous(breaks = unique(summary_long_s4$p)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "VSAT Performance by Item Population Size",
    subtitle = expression("n = 36, set size = 10," ~ rho[1] == 0.7 ~
                            "," ~ rho[2] == 0 ~ "," ~ rho[3] == 0),
    x        = "Item Population Size (p)",
    y        = NULL,
    color    = "Method"
  ) +
  vsat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------------------------------------
# Sample Size
# -------------------------------------------------------

summary_long_s5 <- make_summary(readRDS("results_sample_size.rds"), n)

ggplot(summary_long_s5, aes(x = n, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, ncol = 2) +
  vsat_colors +
  scale_x_continuous(breaks = unique(summary_long_s5$n)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "VSAT Performance by Sample Size",
    subtitle = expression("p = 100, set size = 10," ~ rho[1] == 0.7 ~
                            "," ~ rho[2] == 0 ~ "," ~ rho[3] == 0),
    x        = "Samples per Group (n)",
    y        = NULL,
    color    = "Method"
  ) +
  vsat_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------------------------------------
# Null Simulation - FDR Control
# -------------------------------------------------------

results_null <- readRDS("results_null.rds")

summary_null <- results_null |>
  group_by(rho_1, method) |>
  summarise(
    mean_false_discoveries = mean(returned_size, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(method = case_when(
    method == "kw"      ~ "Kruskal-Wallis",
    method == "permute" ~ "Permutation"
  ))

ggplot(summary_null, aes(x = rho_1, y = mean_false_discoveries,
                         color = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 0.5) +
  vsat_colors +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.1),
    labels = seq(0, 0.9, by = 0.1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "False Discovery Counts Under Non-Differential Correlation",
    subtitle = "p = 100, n = 36, set size = 10",
    y        = "Mean size of incorrect set",
    x        = "Background Correlation Strength",
    color    = "Method"
  ) +
  vsat_theme()
