# ------------------------------------------------------ #
# ------------------- Testing Function ----------------- #
# ------------------------------------------------------ #

source("helpers.R")
library(tidyverse)

VSAT <- function(data,
                 A0      = character(0),
                 method  = c("kw", "permute"),
                 alpha   = 0.05,
                 n_perm  = 1000,
                 K       = 5,
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
    
    pval_fn <- function(A) compute_kw_pvalues(data, uvw, A, taxa)
    
    if (length(A0) == 0) {
      A <- initialize_A(score_fn, pval_fn, taxa, K = K, M = M)
    } else {
      A <- A0
    }
  }
  
  if (method == "permute") {
    
    # Observed correlation matrices
    corr_matrices <- lapply(data, safe_cor)
    taxa          <- Reduce(intersect, lapply(corr_matrices, rownames))
    corr_matrices <- lapply(corr_matrices, function(m) m[taxa, taxa])
    
    # --- KEY CHANGE ---
    # Store only column-group assignments (tiny) instead of 600 pre-computed
    # correlation matrices (≈1.2 GB at p = 500).  .compute_pvals_lean()
    # rebuilds the three group matrices on the fly for each permutation.
    all_cols       <- do.call(cbind, data)
    grp_labels     <- rep(seq_along(data), sapply(data, ncol))
    perm_assignments <- replicate(n_perm, sample(grp_labels), simplify = FALSE)
    
    score_fn <- function(j) {
      means <- sapply(corr_matrices, function(mat) {
        mean(mat[j, colnames(mat) != j])
      })
      var(means)
    }
    
    pval_fn <- function(A) {
      .compute_pvals_lean(corr_matrices, all_cols, perm_assignments,
                          A, taxa, n_perm)
    }
    
    if (length(A0) == 0) {
      scores <- sapply(taxa, score_fn)
      A      <- names(sort(scores, decreasing = TRUE))[1:min(K, length(taxa))]
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