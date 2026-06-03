# ---------------------------------------------- #
# ------------- Testing Generation ------------- #
# ---------------------------------------------- #

library(tidyverse)
library(MASS)
library(mvtnorm)

treatment_data_generation_3d <- function(rows, cols, A_cov, B_cov, C_cov, test_group) {
  
  make_sigma <- function(cov_val) {
    sigma <- diag(rows)
    sigma[test_group, test_group] <- cov_val
    diag(sigma) <- 1
    sigma
  }
  
  lapply(c(A_cov, B_cov, C_cov), function(cov_val) {
    samples <- t(mvrnorm(cols, rep(0, rows), make_sigma(cov_val)))
    data.frame(samples, row.names = paste0("V", 1:rows))
  })
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
