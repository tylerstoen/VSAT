library(tidyverse)
library(MASS)

class_dat <- read.csv(here::here("vineyard/rarified_counts_by_Class.csv"))

class_dat <- class_dat %>% 
  mutate(
    sample = paste(Year, Season, Replicate)
  )


# ------------------------------------------------
# Crop Cover
# ------------------------------------------------

cc_pivoted <- class_dat %>% 
  filter(Experiment == "Crop Cover") %>% 
  group_by(Class, sample, Treatment) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = total,
    values_fill = 0
  )

cc_dfs <- split(cc_pivoted, cc_pivoted$Treatment)

cc_matrices <- lapply(cc_dfs, function(df) {
  df %>%
    column_to_rownames("Class") %>%
    .[, -1]                           
})

cc_matrices_clean <- lapply(cc_matrices, function(mat) {
  keep <- apply(mat, 1, function(x) sd(x) > 0 & !all(is.na(x)))
  mat[keep, , drop = FALSE]
})

# keep only taxa present in all groups
common_taxa <- Reduce(intersect, lapply(cc_matrices_clean, rownames))
cc_matrices_clean <- lapply(cc_matrices_clean, function(mat) {
  mat[common_taxa, , drop = FALSE]
})

lapply(cc_matrices_clean, ncol)

cc_k.w_test <- VSAT(cc_matrices_clean, method = "kw")
cc_perm_test <- VSAT(cc_matrices_clean, method = "permute", n_perm=200)

# ---------------------------------------------
# Fertilizer
# ---------------------------------------------

fert_pivoted <- class_dat %>% 
  filter(Experiment == "Fertilizer") %>% 
  group_by(Class, sample, Treatment) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = total,
    values_fill = 0
  )

fert_dfs <- split(fert_pivoted, fert_pivoted$Treatment)

fert_matrices <- lapply(fert_dfs, function(df) {
  df %>%
    column_to_rownames("Class") %>%
    .[, -1]                           
})

zero_var_rows <- unique(unlist(lapply(fert_matrices, function(mat) {
  rownames(mat)[apply(mat, 1, var) == 0]
})))

fert_matrices_clean <- lapply(fert_matrices, function(mat) {
  mat[!(rownames(mat) %in% zero_var_rows), , drop = FALSE]
})


fert_k.w_test <- VSAT(fert_matrices_clean, method = "kw")
fert_perm_test <- VSAT(fert_matrices_clean, method = "permute", n_perm=200)


# -----------------------------------------------------------
# Capstone Modules Starting Points
# -----------------------------------------------------------

run_vsat_modules <- function(data, modules, method = "kw",
                             n_perm = 200, alpha = 0.05) {
  
  results <- lapply(modules, function(module) {
    returned_set <- VSAT(data, A0 = module, method = method,
                         n_perm = n_perm, alpha = alpha)
    
    data.frame(
      starting_set  = paste(sort(module),       collapse = ", "),
      returned_set  = paste(sort(returned_set), collapse = ", "),
      returned_size = length(returned_set)
    )
  })
  
  bind_rows(results) |>
    group_by(returned_set, returned_size) |>
    summarise(
      starting_sets = paste(starting_set, collapse = " | "),
      n_modules     = n(),
      .groups       = "drop"
    ) |>
    arrange(desc(returned_size))
}


# -------------
# Crop Cover
# -------------

HW_LW <- readRDS("capstone_modules/HW-LW.rds")
HW_NCC <- readRDS("capstone_modules/HW-NCC.rds")
LW_HW <- readRDS("capstone_modules/LW-HW.rds")
LW_NCC <- readRDS("capstone_modules/LW-NCC.rds")
NCC_HW <- readRDS("capstone_modules/NCC-HW.rds")
NCC_LW <- readRDS("capstone_modules/NCC-LW.rds")


# combine all pairwise modules into one list
all_modules_cc <- c(HW_LW, HW_NCC, LW_HW, LW_NCC, NCC_HW, NCC_LW)

# run VSAT on all of them
results_cc <- run_vsat_modules(
  data   =  cc_matrices_clean,
  modules = all_modules,
  method = "permute"
)

saveRDS(results_cc, "capstone_mods_to_VSAT.rds")


# ---------------
# Fertilizer
# ---------------

NF_OF <- readRDS("capstone_modules/NF-OF.rds")
NF_SF <- readRDS("capstone_modules/NF-SF.rds")
OF_NF <- readRDS("capstone_modules/OF-NF.rds")
OF_SF <- readRDS("capstone_modules/OF-SF.rds")
SF_NF <- readRDS("capstone_modules/SF-NF.rds")
SF_OF <- readRDS("capstone_modules/SF-OF.rds")


# combine all pairwise modules into one list
all_modules_fert <- c(NF_OF, NF_SF, OF_NF, OF_SF, SF_NF, SF_OF)

# run VSAT on all of them
results_fert <- run_vsat_modules(
  data   =  fert_matrices_clean,
  modules = all_modules_fert,
  method = "permute"
)

saveRDS(results_cc, "capstone_mods_to_VSAT.rds")


# ---------------------------------------------
# Heatmaps of Result
# ---------------------------------------------

library(gt)

targets <- c("caldilineae", "sphingobacteriia")

cor_mats <- lapply(fert_matrices_clean, safe_cor)
names(cor_mats) <- c("NF", "OF", "SF")

sapply(cor_mats, function(mat) mat["caldilineae", "sphingobacteriia"]) |>
  as.data.frame() |>
  rownames_to_column("Treatment") |>
  rename(Correlation = 2) |>
  gt() |>
  fmt_number(columns = Correlation, decimals = 3) |>
  cols_label(
    Treatment   = "Treatment",
    Correlation = "Correlation (Caldilineae, Sphingobacteriia)"
  ) |>
  tab_header(
    title = "Pairwise Correlation Across Fertilizer Treatments"
  )
