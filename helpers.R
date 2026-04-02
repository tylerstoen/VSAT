# ---------------------------------------------------- #
# ------------------- Helpers ------------------------ #
# ---------------------------------------------------- #

# --------------------
# General Computations
# --------------------

compute_S <- function(corr_matrix, j, A) { # where j is the target and A is a set
  A <- setdiff(A, j)
  if (length(A) == 0) return(0)
  return(1/length(A) * (sum(corr_matrix[j,A])))
}

compute_delta <- function(corr_matrices, j, A) {
  s_values <- sapply(corr_matrices, 
                     function(corr_matrix) compute_S(corr_matrix, j, A))
  avg_s <- 1/length(corr_matrices) * sum(s_values)
  return(sum((s_values-avg_s)^2))
}

# compute m vectors
compute_m_vecs <- function(uvw_vals, A, j) {
  
  n_treat <- length(uvw_vals)
  
  # helper: compute column means safely
  col_mean_or_zero <- function(mat, rows) {
    if (length(rows) == 0) {
      rep(0, ncol(mat))
    } else {
      colMeans(mat[rows, , drop = FALSE])
    }
  }
  
  if (j %in% A) {
    A_use <- setdiff(A, j)
  } else {
    A_use <- A
  }
  
  m1 <- col_mean_or_zero(uvw_vals[[1]], A_use)
  m2 <- col_mean_or_zero(uvw_vals[[2]], A_use)
  m3 <- col_mean_or_zero(uvw_vals[[3]], A_use)
  
  list(
    m1 = m1,
    m2 = m2,
    m3 = m3
  )
}


compute_uvw <- function(data) {
  uvw_vals <- lapply(data, function(mat){
    t(scale(t(mat)))
  })
}






# ---------------------------------------------------------------------------- #



# ----------------------------
# Permutation Specific Helpers
# ----------------------------




get_permutations <- function(data) {
  n_conditions <- length(data)
  n_row <- nrow(data[[1]])
  n_samples <- ncol(data[[1]][-c(1,2)])
  classes <- data[[1]]$Class
  sample_names <- colnames(data[[1]])[-c(1,2)]
  
  # build permuted dataframes
  permuted_data <- vector("list", n_conditions)
  for (i in seq_len(n_conditions)) {
    permuted_data[[i]] <- data[[i]]
  }
  
  # loop over variables and samples
  for (row in seq_len(n_row)) {
    for (sample_j in seq_len(n_samples)) {
      # extract counts across conditions for this (var, sample)
      counts <- sapply(data, function(df) df[row, sample_names[sample_j]])
      
      # permute across conditions
      permuted_counts <- sample(counts, length(counts))
      
      # save into permuted dataframes
      for (cond_k in seq_len(n_conditions)) {
        permuted_data[[cond_k]][row, sample_names[sample_j]] <- permuted_counts[cond_k]
      }
    }
  }
  return(permuted_data)
}


permute_data <- function(data, nperm) {
  permed_corr_mats <- vector("list", nperm)
  
  for (n in 1:nperm) {
    orig_order <- c(rep("A", ncol(data[[1]])), rep("B", ncol(data[[2]])), rep("C", ncol(data[[3]])))
    all_data_together <- cbind(data[[1]], data[[2]], data[[3]])
    
    perm_1 <- sample(orig_order)
    
    fake_A <- all_data_together[, perm_1 == "A"]
    fake_B <- all_data_together[, perm_1 == "B"]
    fake_C <- all_data_together[, perm_1 == "C"]
    
    rownames(fake_A) <- rownames(data[[1]])
    rownames(fake_B) <- rownames(data[[2]])
    rownames(fake_C) <- rownames(data[[3]])
    
    permed_data <- list(fake_A, fake_B, fake_C)
    permed_corr_mats[[n]] <- lapply(permed_data, function(mat) cor(t(mat)))
  }
  return(permed_corr_mats)
} # returns a list of three correlation matrices





vineyard_data_to_matrix <- function(data) { #input list of data (object from split() works)
  names(data) <- seq_along(data)
  
  # Three dataframe inputs w/ variable, group, and all samples as columns
  mats <- lapply(data, function(df){
    df %>%
      column_to_rownames(var = colnames(df)[1]) %>%  # first column as rownames
      .[, -1]
  }) # drop variable and group cols
  
  zero_var_rows <- unique(unlist(lapply(mats, function(mat) {
    rownames(mat)[apply(mat, 1, var) == 0]
  })))
  
  clean_mats <- lapply(mats, function(mat) {
    mat[!(rownames(mat) %in% zero_var_rows), , drop = FALSE]
  })
  
  return(clean_mats)
}





