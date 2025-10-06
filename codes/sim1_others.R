suppressPackageStartupMessages(library(np))
suppressPackageStartupMessages(library(vcrpart))
library(ranger)
library(BART)
source("common_sim.R")
source("competitor_wrappers.R")
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
stopifnot(!is.na(seed), seed >= 0)

######### (1)  KS #########
run_ks <- function(seed = 1,
                   n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                   p = 3, R = 20, sigma = 1.0, B = 50) {
  
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  res <- kernel_smoothing_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z, B = B)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  a_mse <- NA; b_cov <- NA
  if (!is.null(res$test$beta)) {
    bt <- res$test$beta
    dn <- dimnames(bt)[[2]]
    if (!is.null(dn) && all(c("MEAN","L95","U95") %in% dn)) {
      mean_hat <- as.numeric(bt[, "MEAN", , drop = FALSE])
      L        <- as.numeric(bt[, "L95",  , drop = FALSE])
      U        <- as.numeric(bt[, "U95",  , drop = FALSE])
      truth_at <- as.numeric(sapply(truth_list, function(f) f(te$Z[, 1:5, drop = FALSE])))
      
      a_mse <- mean((mean_hat - truth_at)^2)
      b_cov <- mean((truth_at >= L) & (truth_at <= U))
    }
  }
  
  list(a_mse_beta = a_mse,
       b_cov_beta = b_cov,
       c_rmse_pred = unname(pred["rmse"]),
       d_cov_pred  = unname(pred["cov"]))
}

res_KS <- run_ks(
  seed = seed,
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0, B = 50
)

######### (2)  BTVCM #########
run_btvcm <- function(seed = 1,
                      n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                      p = 3, R = 20, sigma = 1.0,
                      lambda = 0.05, ntree = 200, B = 0) {
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  res <- boosted_tvcm_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z,
                              intercept = TRUE, lambda = lambda, ntree = ntree, B = B)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  a_mse <- NA; b_cov <- NA
  if (!is.null(res$test$beta)) {
    bt <- res$test$beta
    dn <- dimnames(bt)[[2]]
    if (!is.null(dn) && all(c("MEAN","L95","U95") %in% dn)) {
      mean_hat <- bt[, "MEAN", , drop = FALSE]
      L        <- bt[, "L95",  , drop = FALSE]
      U        <- bt[, "U95",  , drop = FALSE]
      truth_at <- sapply(truth_list, function(f) f(te$Z[, 1:5, drop = FALSE]))
      mean_hat <- as.numeric(mean_hat)
      L        <- as.numeric(L)
      U        <- as.numeric(U)
      truth_at <- as.numeric(truth_at)
      a_mse <- mean((mean_hat - truth_at)^2)
      b_cov <- mean((truth_at >= L) & (truth_at <= U))
    }
  }
  
  list(a_mse_beta = a_mse,
       b_cov_beta = b_cov,
       c_rmse_pred = unname(pred["rmse"]),
       d_cov_pred  = unname(pred["cov"]))
}

res_BTVCM <- run_btvcm(
  seed = seed,
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0, B = 50
)
                         
######### (3) TVC #########
run_tvc <- function(seed = 1,
                    n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                    p = 3, R = 20, sigma = 1.0,
                    B = 0) {  # B = bootstrap reps for beta intervals
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  res  <- tvc_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z, B = B)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  a_mse <- NA; b_cov <- NA
  if (!is.null(res$test$beta)) {
    bt <- res$test$beta  # [n_test x {MEAN, L95, U95, ...} x (p+1)]
    dn <- dimnames(bt)[[2]]
    if (!is.null(dn) && all(c("MEAN", "L95", "U95") %in% dn)) {
      mean_hat <- as.numeric(bt[, "MEAN", , drop = FALSE])
      L        <- as.numeric(bt[, "L95",  , drop = FALSE])
      U        <- as.numeric(bt[, "U95",  , drop = FALSE])
      truth_at <- as.numeric(sapply(truth_list, function(f) f(te$Z[, 1:5, drop = FALSE])))
      
      a_mse <- mean((mean_hat - truth_at)^2)
      b_cov <- mean((truth_at >= L) & (truth_at <= U))
    }
  }
  
  list(a_mse_beta = a_mse,
       b_cov_beta = b_cov,
       c_rmse_pred = unname(pred["rmse"]),
       d_cov_pred  = unname(pred["cov"]))
}

res_TVC <- run_tvc(
  seed = seed,
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0, B = 50
)

######### (4) LM #########
              
run_lm <- function(seed = 1,
                   n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                   p = 3, R = 20, sigma = 1.0
) {
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  res <- lm_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  if (!is.null(res$test$beta)) {
    bt <- res$test$beta
    dn <- dimnames(bt)[[2]]
    if (all(c("MEAN","L95","U95") %in% dn)) {
      mean_hat <- bt[, "MEAN", ]
      truth_at <- sapply(truth_list, function(f) f(te$Z[,1:5, drop=FALSE]))
      a_mse <- mean((mean_hat - truth_at)^2)
      # coverage via provided L95/U95
      L <- bt[, "L95", ]; U <- bt[, "U95", ]
      b_cov <- mean((truth_at >= L) & (truth_at <= U))
    } else { a_mse <- NA; b_cov <- NA }
  } else { a_mse <- NA; b_cov <- NA }
  
  list(a_mse_beta = a_mse, b_cov_beta = b_cov, c_rmse_pred = pred["rmse"], d_cov_pred = pred["cov"])
}

res_LM <- run_lm(
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0,seed=seed
)                     

######### (5) BART #########
run_bart <- function(seed = 1,
                        n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                        p = 3, R = 20, sigma = 1.0,
                        nd = 2000, burn = 400
) {
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  N_train <<- nrow(tr$X); N_test <<- nrow(te$X); Y_test <<- te$Y
  
  res <- bart_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z, nd = nd, burn = burn, thin = 1)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  list(a_mse_beta = NA_real_, b_cov_beta = NA_real_, c_rmse_pred = pred["rmse"], d_cov_pred = pred["cov"])
}

res_BART <- run_bart(
    n_subj_train = 1000, n_subj_test = 200, m_per = 1,
    p = 3, R = 20, sigma = 1.0, seed=seed
)

######### (6) ERT #########
run_ert <- function(seed = 1,
                       n_subj_train = 1000, n_subj_test = 250, m_per = 1,
                       p = 3, R = 20, sigma = 1.0
) {
  set.seed(seed)
  tr <- gen_data(n_subj_train, m_per, p, R, sigma)
  te <- gen_data(n_subj_test,  m_per, p, R, sigma)
  
  N_train <<- nrow(tr$X); N_test <<- nrow(te$X); Y_test <<- te$Y
  
  res <- extraTrees_wrapper(tr$Y, tr$X, tr$Z, te$X, te$Z)
  pred <- pred_from_wrapper(res$test$ystar, te$Y)
  
  list(a_mse_beta = NA_real_, b_cov_beta = NA_real_, c_rmse_pred = pred["rmse"], d_cov_pred = pred["cov"])
}

res_ERT <- run_ert(
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0, seed=seed
)                            
