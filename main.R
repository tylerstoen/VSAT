
# ------------------------------------------------------ #
# ------------------- Testing Function ----------------- #
# ------------------------------------------------------ #
source("helpers.R")
library(tidyverse)

A_star_search <- function(data,
                          A0      = character(0),
                          method  = c("kw", "permute"),
                          alpha   = 0.05,
                          n_perm  = 1000,
                          K       = 10,
                          M       = 20) {
  
  method <- match.arg(method)
  reps   <- 0
  A_history <- list()
  
  if (method == "kw") {
    
    taxa <- rownames(data[[1]])
    uvw  <- compute_uvw(data)
    n1   <- ncol(data[[1]]); n2 <- ncol(data[[2]]); n3 <- ncol(data[[3]])
    
    score_fn <- function(j) {
      means <- sapply(uvw, function(mat) mean(mat[j, ]))
      var(means)
    }
    
    pval_fn <- function(A) compute_kw_pvals(uvw, A, taxa, n1, n2, n3)
    
    if (length(A0) == 0) {
      A <- initialize_A(score_fn, pval_fn, taxa, K = K, M = M)
    } else {
      A <- A0
    }
  }
  
  if (method == "permute") {
    
    corr_matrices        <- lapply(data, function(mat) cor(t(mat)))
    taxa                 <- rownames(corr_matrices[[1]])
    permed_corr_matrices <- permute_data(data, n_perm)
    
    score_fn <- function(j) {
      means <- sapply(corr_matrices, function(mat) {
        mean(mat[j, colnames(mat) != j])
      })
      var(means)
    }
    
    pval_fn <- function(A) {
      compute_permuted_pvalues(corr_matrices, permed_corr_matrices, A, taxa, n_perm)
    }
    
    if (length(A0) == 0) {
      scores <- sapply(taxa, score_fn)
      A <- names(sort(scores, decreasing = TRUE))[1:min(K, length(taxa))]
    } else {
      A <- A0
    }
  }
  
  repeat {
    
    p_values <- pval_fn(A)
    
    p.adj <- p.adjust(p_values, method = "BH")
    sig   <- names(p_values)[!is.na(p.adj) & p.adj < alpha]
    
    if (length(A) == 1) sig <- union(sig, A)
    if (setequal(A, sig)) break
    
    sig_str         <- paste(sort(sig), collapse = ",")
    history_strings <- sapply(A_history, paste, collapse = ",")
    
    if (sig_str %in% history_strings) {
      message("Exiting early: detected a repeating pattern in A.")
      break
    } else {
      A_history[[length(A_history) + 1]] <- sort(sig)
    }
    
    print(A)
    A    <- sig
    reps <- reps + 1
    print(reps)
  }
  
  return(A)
}


# -------------------------------------------------------- #
# ----------------- K-W test computation ----------------- #
# -------------------------------------------------------- #

compute_kw_pvalues <- function(data, uvw, A, taxa) {
  
  n1 <- length(data[[1]])
  n2 <- length(data[[2]])
  n3 <- length(data[[3]])
  
  pvals <- numeric(length(taxa))
  names(pvals) <- taxa
  
  for (j in taxa) {
    
    M <- compute_m_vecs(uvw, A, j)
    
    m1 <- M$m1
    m2 <- M$m2
    m3 <- M$m3
    
    u_tild <- uvw[[1]][j,] * m1
    v_tild <- uvw[[2]][j,] * m2
    w_tild <- uvw[[3]][j,] * m3
    
    x <- c(u_tild, v_tild, w_tild)
    g <- c(rep(1, n1), rep(2, n2), rep(3, n3))
    
    pvals[j] <- kruskal.test(x = x, g = g)$p.value
  }
  
  return(pvals)
}

# ----------------------------------------------------- #
# --------------- p-vals from permuation -------------- #
# ----------------------------------------------------- #

compute_permuted_pvalues <- function(corr_matrices, permed_corr_mats_list, A, taxa, n_perm) {
  
  delta_values <- numeric(length(taxa))
  names(delta_values) <- taxa
  
  for (tax in taxa) {
    if (tax %in% A) {
      delta_values[tax] <- compute_delta(corr_matrices, tax, setdiff(A, tax))
    } else {
      delta_values[tax] <- compute_delta(corr_matrices, tax, A)
    }
  }
  
  delta_perm <- matrix(NA, nrow = length(taxa), ncol = n_perm)
  rownames(delta_perm) <- taxa
  
  for (perm in seq_len(n_perm)) {
    perm_mats <- permed_corr_mats_list[[perm]]
    
    # skip permutations where any group lost taxa due to zero SD
    available_taxa <- Reduce(intersect, lapply(perm_mats, rownames))
    
    delta_perm[available_taxa, perm] <- sapply(available_taxa, function(j) {
      compute_delta(perm_mats, j, A)
    })
  }
  
  # NA columns (bad permutations) are ignored by mean()
  p_values <- sapply(seq_along(delta_values), function(i) {
    mean(delta_perm[i, ] >= delta_values[i], na.rm = TRUE)
  })
  
  names(p_values) <- taxa
  return(p_values)
}

# --------------------------------------------------------- #
# ----------- Generalized Initialization Procedure -------- #
# --------------------------------------------------------- #

# score_fn : taxa -> numeric scalar  (higher = more promising start)
# pval_fn  : A    -> named p-value vector
initialize_A <- function(score_fn, pval_fn, taxa, K = 10, M = 20) {
  
  scores        <- sapply(taxa, score_fn)
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



