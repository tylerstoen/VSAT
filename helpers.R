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

# -------------------------------------------------------
# Vectorised delta for ALL taxa at once 
# -------------------------------------------------------
.deltas_fast <- function(mats, taxa, A) {
  A_use <- intersect(A, taxa)          # anchors actually present
  n_g   <- length(mats)
  
  if (length(A_use) == 0) return(setNames(rep(0, length(taxa)), taxa))
  
  # S_mat[i, g] = mean correlation of taxa[i] with anchor set A (self-excluded)
  S_mat <- matrix(0, nrow = length(taxa), ncol = n_g,
                  dimnames = list(taxa, NULL))
  
  for (g in seq_len(n_g)) {
    cm         <- mats[[g]]
    S_mat[, g] <- rowMeans(cm[taxa, A_use, drop = FALSE])  # bulk vectorised op
  }
  
  # Correct the ≤ K taxa that appear in A themselves (self-exclusion)
  for (j in intersect(taxa, A_use)) {
    A_no_j <- setdiff(A_use, j)
    for (g in seq_len(n_g))
      S_mat[j, g] <- if (length(A_no_j)) mean(mats[[g]][j, A_no_j]) else 0
  }
  
  avg <- rowMeans(S_mat)
  rowSums((S_mat - avg)^2)             # sum of squared deviations across groups
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


safe_cor <- function(mat) {
  keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  if (sum(keep) < 2) return(NULL)
  cor(t(mat[keep, , drop = FALSE]))
}


# ----------------------------------------------------- #
# --------------- p-vals from permutation ------------- #
# ----------------------------------------------------- #

compute_permuted_pvalues <- function(corr_matrices, permed_corr_mats_list,
                                     A, taxa, n_perm) {
  
  delta_obs  <- .deltas_fast(corr_matrices, taxa, A)
  delta_perm <- matrix(NA_real_, nrow = length(taxa), ncol = n_perm,
                       dimnames = list(taxa, NULL))
  
  for (perm in seq_len(n_perm)) {
    pm    <- permed_corr_mats_list[[perm]]
    avail <- Reduce(intersect, lapply(pm, rownames))
    ct    <- intersect(taxa, avail)
    if (length(ct) < 2) next
    pm_sub           <- lapply(pm, function(m) m[ct, ct, drop = FALSE])
    delta_perm[ct, perm] <- .deltas_fast(pm_sub, ct, A)
  }
  
  p_vals <- rowMeans(sweep(delta_perm, 1, delta_obs, ">="), na.rm = TRUE)
  names(p_vals) <- taxa
  p_vals
}


# --------------------------------------------------------------------- #
# Memory-efficient permuted p-values (used when method == "permute")    #
# Stores only column-group assignments instead of 600 corr matrices,    #
# keeping RAM near-constant regardless of p.                             #
# Called from inside VSAT; not a drop-in for compute_permuted_pvalues.  #
# --------------------------------------------------------------------- #
.compute_pvals_lean <- function(corr_matrices, all_cols, perm_assignments,
                                A, taxa, n_perm) {
  
  delta_obs    <- .deltas_fast(corr_matrices, taxa, A)
  exceed_count <- setNames(integer(length(taxa)), taxa)
  valid_count  <- setNames(integer(length(taxa)), taxa)
  
  for (perm in seq_len(n_perm)) {
    grp <- perm_assignments[[perm]]
    pm  <- lapply(sort(unique(grp)),
                  function(g) safe_cor(all_cols[, grp == g, drop = FALSE]))
    
    if (any(sapply(pm, is.null))) next
    avail <- Reduce(intersect, lapply(pm, rownames))
    ct    <- intersect(taxa, avail)
    if (length(ct) < 2) next
    
    pm_sub <- lapply(pm, function(m) m[ct, ct, drop = FALSE])
    d_perm <- .deltas_fast(pm_sub, ct, A)
    
    idx               <- match(ct, taxa)
    exceed_count[idx] <- exceed_count[idx] + (d_perm >= delta_obs[idx])
    valid_count[idx]  <- valid_count[idx] + 1L
    
    rm(pm, pm_sub)
  }
  
  p_vals <- exceed_count / pmax(valid_count, 1L)
  names(p_vals) <- taxa
  p_vals
}


# -------------------------------------
# Generalized Initialization Procedure 
# -------------------------------------

initialize_A <- function(score_fn, pval_fn, taxa, K = 10, M = 20) {
  
  scores         <- sapply(taxa, score_fn)
  top_candidates <- names(sort(scores, decreasing = TRUE))[1:min(M, length(taxa))]
  
  best_A     <- NULL
  best_score <- Inf
  
  for (start in top_candidates) {
    
    pvals       <- pval_fn(start)
    A_candidate <- names(sort(pvals))[1:min(K, length(pvals))]
    
    cand_pvals  <- pval_fn(A_candidate)
    score       <- mean(sort(cand_pvals)[1:min(K, length(cand_pvals))])
    
    if (score < best_score) {
      best_score <- score
      best_A     <- A_candidate
    }
  }
  
  return(best_A)
}


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