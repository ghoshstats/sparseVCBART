library(MASS)
library(sparseVCBART)
library(VCBART)
source("common_sim.R")

set.seed(1)
run_sim <- function(
    n_subj_train = 1000, n_subj_test = 250, m_per = 1,
    p = 3, R = 20, sigma = 1.0,
    M_sparse = 50, M_van = 50, nd = 1200, burn = 400,
    n_cuts = 50, grid_g = 200,
    slope_tau_mult = NULL    
) {
  ## ----- generate ONE training set and ONE test set  -----
  train <- gen_data(n_subj = n_subj_train, m_per = m_per, p = p, R = R, sigma = sigma)
  test  <- gen_data(n_subj = n_subj_test,  m_per = m_per, p = p, R = R, sigma = sigma)
  
  Zc_tr <- train$Z; Zk_tr <- matrix(0L, nrow = 1, ncol = 1)
  Zc_te <- test$Z;  Zk_te <- matrix(0L, nrow = 1, ncol = 1)
  
  # cutpoints from training modifiers only
  cp_tr <- make_cutpoints(Zc_tr, n_cuts = n_cuts)
  
  ## ----- sparseVCBART ------
  if (!is.null(slope_tau_mult)) {
    tau0  <- 0.5 / sqrt(M_sparse)
    tau_s <- c(tau0, rep(slope_tau_mult * tau0, p))
  } else {
    tau_s <- NULL
  }
  fit_sparse <- sparseVCBART::VCBART_ind(
    Y_train = train$Y,
    subj_id_train = train$subj_id,
    ni_train = train$ni,
    X_train = train$X,
    Z_cont_train = Zc_tr,
    Z_cat_train = Zk_tr,
    cutpoints_list = cp_tr,
    M = M_sparse, nd = nd, burn = burn,
    tau = tau_s,
    save_samples = TRUE, save_trees = TRUE, verbose = FALSE
  )
  
  ## ----- vanilla VCBART  -----
  tau_vec <- rep(0.5 / sqrt(M_van), p + 1)
  fit_van <- VCBART::VCBART_ind(
    Y_train = train$Y,
    subj_id_train = train$subj_id,
    ni_train = train$ni,
    X_train = train$X,
    Z_cont_train = Zc_tr,
    Z_cat_train = Zk_tr,
    cutpoints_list = cp_tr,
    tau = tau_vec,
    M = M_van, nd = nd, burn = burn,
    save_samples = TRUE, save_trees = TRUE, verbose = FALSE
  )
  
  Zg <- matrix(0.5, nrow = grid_g, ncol = R)
  Zg[, 1] <- seq(0, 1, length.out = grid_g)
  Zg[, 2] <- seq(0, 1, length.out = grid_g)
  Zk_g <- matrix(0L, nrow = grid_g, ncol = 1)
  
  b_sparse_g <- sparseVCBART::predict_betas(fit_sparse, Z_cont = Zg, Z_cat = Zk_g, verbose = FALSE)
  b_van_g    <- VCBART::predict_betas(fit_van,    Z_cont = Zg, Z_cat = Zk_g, verbose = FALSE)
  
  sum_sparse <- sparseVCBART::summarize_beta(b_sparse_g)   # [g x {MEAN,L95,U95} x (p+1)]
  sum_van    <- VCBART::summarize_beta(b_van_g)
  
  # (a) Average MSE #
  mse_s <- mse_beta_on_grid(sum_sparse, Zg, truth_list)
  mse_v <- mse_beta_on_grid(sum_van,    Zg, truth_list)
  
  # (b) Average 95% coverage #
  cov_s <- cov_beta_on_grid(sum_sparse, Zg, truth_list, level = 0.95)
  cov_v <- cov_beta_on_grid(sum_van,    Zg, truth_list, level = 0.95)
  
  ## ----- (c,d) predictive metrics on the test set -----
  b_sparse_te <- sparseVCBART::predict_betas(fit_sparse, Z_cont = Zc_te, Z_cat = Zk_te, verbose = FALSE)
  b_van_te    <- VCBART::predict_betas(fit_van,    Z_cont = Zc_te, Z_cat = Zk_te, verbose = FALSE)
  
  pred_s <- pred_metrics(b_sparse_te, fit_sparse$sigma, test$X, test$Y, level = 0.95)
  pred_v <- pred_metrics(b_van_te,    fit_van$sigma,    test$X, test$Y, level = 0.95)
  
  list(
    a_mse_beta = c(sparse = mse_s$avg, vanilla = mse_v$avg),
    b_cov_beta = c(sparse = cov_s$avg, vanilla = cov_v$avg),
    c_rmse_pred = c(sparse = pred_s$rmse, vanilla = pred_v$rmse),
    d_cov_pred  = c(sparse = pred_s$cov,  vanilla = pred_v$cov),
    per_j = list(
      mse_sparse = mse_s$per_j, mse_van = mse_v$per_j,
      cov_sparse = cov_s$per_j, cov_van = cov_v$per_j
    )
  )
}

############ Running a single chain ###########
for(i in 1:25){
  set.seed(i)
  res <- run_sim(
    n_subj_train = 1000, n_subj_test = 200, m_per = 1,
    p = 3, R = 20, sigma = 1.0,
    M_sparse = 50, M_van = 50, nd = 2000, burn = 400,
    n_cuts = 50, grid_g = 200,
    slope_tau_mult = NULL  # or try 2.0 if you want looser sparse slopes
  )
  
  
  print(res$a_mse_beta)     # (a) Average MSE for beta(z)
  print(res$b_cov_beta)     # (b) Average 95% coverage for beta(z)
  print(res$c_rmse_pred)    # (c) Predictive RMSE
  print(res$d_cov_pred)     # (d) 95% prediction interval coverage
}

############# Running multiple chains ################
######################################################
run_sim_chain <- function(chain_seed, ...) {
  set.seed(chain_seed)
  res <- run_sim(...)
  data.frame(
    chain_seed     = chain_seed,
    a_mse_sparse   = as.numeric(res$a_mse_beta["sparse"]),
    a_mse_van      = as.numeric(res$a_mse_beta["vanilla"]),
    b_cov_sparse   = as.numeric(res$b_cov_beta["sparse"]),
    b_cov_van      = as.numeric(res$b_cov_beta["vanilla"]),
    c_rmse_sparse  = as.numeric(res$c_rmse_pred["sparse"]),
    c_rmse_van     = as.numeric(res$c_rmse_pred["vanilla"]),
    d_cov_sparse   = as.numeric(res$d_cov_pred["sparse"]),
    d_cov_van      = as.numeric(res$d_cov_pred["vanilla"]),
    row.names = NULL
  )
}

run_seed_multichain <- function(seed,
                                n_chains = 4,
                                # pass-through args to run_sim():
                                n_subj_train = 1000, n_subj_test = 200, m_per = 1,
                                p = 3, R = 20, sigma = 1.0,
                                M_sparse = 50, M_van = 50, nd = 2000, burn = 400,
                                n_cuts = 50, grid_g = 200,
                                slope_tau_mult = NULL) {

  chain_seeds <- seed + seq_len(n_chains)  
  per_chain <- do.call(rbind, lapply(
    chain_seeds,
    function(cs) run_sim_chain(
      chain_seed = cs,
      n_subj_train = n_subj_train, n_subj_test = n_subj_test, m_per = m_per,
      p = p, R = R, sigma = sigma,
      M_sparse = M_sparse, M_van = M_van,
      nd = nd, burn = burn, n_cuts = n_cuts, grid_g = grid_g,
      slope_tau_mult = slope_tau_mult
    )
  ))

  agg <- as.data.frame(colMeans(per_chain[ , -1, drop = FALSE]))
  names(agg) <- "mean"
  agg$metric <- rownames(agg)
  rownames(agg) <- NULL
  list(seed = seed, per_chain = per_chain, summary = agg)
}

run_all_seeds_multichain <- function(seeds = 1:25, n_chains = 4, ...) {
  out_list <- lapply(seeds, function(s) run_seed_multichain(seed = s, n_chains = n_chains, ...))

  seed_summaries <- do.call(rbind, lapply(out_list, function(x) {
    w <- reshape(
      x$summary,
      timevar = "metric",
      idvar   = NULL,
      direction = "wide"
    )
    w$seed <- x$seed
    w <- w[, c(ncol(w), setdiff(seq_len(ncol(w)-1), integer(0)))]
    names(w) <- sub("^mean\\.", "", names(w))
    w
  }))
  list(
    per_seed = seed_summaries,      # each row = one seed, metrics averaged across chains
    details  = out_list             # includes per-chain tables if you want them
  )
}


res_all <- run_all_seeds_multichain(
  seeds = 1:25,
  n_chains = 4,
  n_subj_train = 1000, n_subj_test = 200, m_per = 1,
  p = 3, R = 20, sigma = 1.0,
  M_sparse = 50, M_van = 50, nd = 2000, burn = 400,
  n_cuts = 50, grid_g = 200,
  slope_tau_mult = NULL
)

head(res_all$per_seed)

