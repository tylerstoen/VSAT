# ------------------------------------------- #
# ---------------- Plotting ----------------- #
# ------------------------------------------- #

library(ggplot2)

# -------------------------------------------------------------------------- #
plot_m_scores <- function(data,
                          A0 = character(0),
                          taxa_subset = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(tibble)
  
  uvw_vals <- compute_uvw(data)
  taxa <- rownames(data[[1]])
  
  if (!is.null(taxa_subset)) {
    taxa <- intersect(taxa, taxa_subset)
  }
  
  out <- vector("list", length(taxa))
  names(out) <- taxa
  
  for (j in taxa) {
    
    M <- compute_m_vecs(uvw_vals, A0, j)
    
    S.j1 <- cov(uvw_vals[[1]][j, ], M$m1)
    S.j2 <- cov(uvw_vals[[2]][j, ], M$m2)
    S.j3 <- cov(uvw_vals[[3]][j, ], M$m3)
    
    out[[j]] <- tibble(
      Taxa = j,
      Treatment = factor(c("T1","T2","T3"),
                         levels = c("T1","T2","T3")),
      m_score = c(S.j1, S.j2, S.j3)
    )
  }
  
  m_df <- bind_rows(out)
  
  ggplot(m_df, aes(Treatment, m_score, color = Treatment)) +
    geom_point() +
    geom_line(aes(group = Taxa)) +
    facet_wrap(~ Taxa) +
    theme_minimal()
}



# ---------------------------------------------------------------------------- #

# ------------------------------------------------------------ #
# Plots individual contribution to correlations                #
# Shows skewness for contributions under treatments w high cov #
# ------------------------------------------------------------ #

plot_contributions <- function(data,
                               A0 = character(0),
                               taxa_subset = NULL,
                               show_means = FALSE) {
  
  library(dplyr)
  library(ggplot2)
  library(tibble)
  
  uvw_vals <- compute_uvw(data)
  taxa <- rownames(data[[1]])
  
  if (!is.null(taxa_subset)) {
    taxa <- intersect(taxa, taxa_subset)
  }
  
  out <- vector("list", length(taxa))
  names(out) <- taxa
  
  for (j in taxa) {
    
    M <- compute_m_vecs(uvw_vals, A0, j)
    
    u_vals <- uvw_vals[[1]][j, ] * M$m1
    v_vals <- uvw_vals[[2]][j, ] * M$m2
    w_vals <- uvw_vals[[3]][j, ] * M$m3
    
    out[[j]] <- tibble(
      Taxa = j,
      Treatment = rep(c("T1","T2","T3"),
                      times = c(length(u_vals),
                                length(v_vals),
                                length(w_vals))),
      Value = c(u_vals, v_vals, w_vals)
    )
  }
  
  contributions <- bind_rows(out)
  
  p <- ggplot(contributions,
              aes(Treatment, Value, color = Treatment)) +
    geom_jitter(width = 0.15,
                height = 0,
                size = 1.8,
                alpha = 0.8) +
    facet_wrap(~ Taxa) +
    theme_minimal()
  
  if (show_means) {
    p <- p +
      stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.5,
        fatten = 1,
        linewidth = 0.8
      )
  }
  
  return(p)
}

# -----------
# ANOVA plot
# -----------

get_contributions <- function(data, A0 = character(0)) {
  
  A <- A0
  taxa <- rownames(data[[1]])
  
  uvw_vals <- lapply(data, function(mat) t(scale(t(mat))))
  
  out <- list()
  
  for (j in taxa) {
    m <- compute_m_vecs(uvw_vals, A, j)
    u_contributions <- uvw_vals[[1]][j, ] * m$m1
    v_contributions <- uvw_vals[[2]][j, ] * m$m2
    w_contributions <- uvw_vals[[3]][j, ] * m$m3
    
    n1 <- length(u_contributions)
    n2 <- length(v_contributions)
    n3 <- length(w_contributions)
    
    out[[j]] <- tibble(
      Taxa = j,
      Treatment = rep(c("T1","T2","T3"), times = c(n1, n2, n3)),
      Value = c(u_contributions, v_contributions, w_contributions)
    )
  }
  bind_rows(out)
}


contributions <- get_contributions(treatment_data_generation_3d(100,20,0.8,0,0,c(1:10)), 
                                   c("V1","V2","V3","V4","V5","V6", "V7", "V8", "V9", "V10"))

subset_taxa <- c("V1", "V20")

contributions |>
  filter(Taxa %in% subset_taxa) |>
  ggplot(aes(Treatment, Value, color = Treatment)) +
  geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.8) +
  labs(y="Contribution",
       title="Contribution breakdown of the Correlation between Two Items") +
  facet_wrap(~Taxa) +
  theme_bw()


