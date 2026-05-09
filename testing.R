# ---------------------------------------------- #
# ------------- Testing Generation ------------- #
# ---------------------------------------------- #

library(tidyverse)
library(MASS)
library(mvtnorm)

treatment_data_generation_3d <- function(rows, cols, A_cov, B_cov, C_cov, test_group) {
  sigma_A <- diag(rows)  # makes a 100x100 identity matrix
  sigma_A[test_group, test_group] <- A_cov #gives the first 10 rows positive correlations
  sigma_A[row(sigma_A) == col(sigma_A)] = 1 # makes sure diagonal is one
  
  sigma_B <- diag(rows)
  sigma_B[test_group,test_group] <- B_cov
  sigma_B[row(sigma_B) == col(sigma_B)] = 1
  
  sigma_C <- diag(rows)
  sigma_C[test_group,test_group] <- C_cov
  sigma_C[row(sigma_C) == col(sigma_C)] = 1
  
  sample_A <- t(mvrnorm(cols, rep(0,rows), sigma_A))
  sample_B <- t(mvrnorm(cols, rep(0,rows), sigma_B))
  sample_C <- t(mvrnorm(cols, rep(0,rows), sigma_C))
  
  test_dfs <- list(
    data.frame(sample_A, row.names = paste0("V", 1:rows)), 
    data.frame(sample_B, row.names = paste0("V", 1:rows)), 
    data.frame(sample_C, row.names = paste0("V", 1:rows))
  )
  return(test_dfs)
}

treatment_data_generation_bg <- function(rows, cols, A_cov, B_cov, C_cov, 
                                         test_group, rho_bg = 0.2) {
  
  make_sigma <- function(cov_val) {
    sigma <- matrix(rho_bg, nrow = rows, ncol = rows)  # global background
    diag(sigma) <- 1
    sigma[test_group, test_group] <- cov_val            # DC block
    diag(sigma) <- 1
    sigma
  }
  
  sigma_A <- make_sigma(A_cov)
  sigma_B <- make_sigma(B_cov)
  sigma_C <- make_sigma(C_cov)
  
  list(
    data.frame(t(mvrnorm(cols, rep(0, rows), sigma_A)), 
               row.names = paste0("V", 1:rows)),
    data.frame(t(mvrnorm(cols, rep(0, rows), sigma_B)),
               row.names = paste0("V", 1:rows)),
    data.frame(t(mvrnorm(cols, rep(0, rows), sigma_C)),
               row.names = paste0("V", 1:rows))
  )
}
