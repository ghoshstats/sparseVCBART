## ===================== Simulation 2 (VCBART variants) ======================
suppressPackageStartupMessages({
  library(MASS)
  library(VCBART)
  library(sparseVCBART)
})

## ----- truths -----
beta0 <- function(z)  3*z[,1] + ((2 - 5*(z[,2] > 0.5)) * sin(pi*z[,1])) - 2*(z[,2] > 0.5)
beta1 <- function(z)  ((3 - 3*z[,2]) * (z[,1] > 0.6) - 10*sqrt(pmax(z[,1], 0)) * (z[,1] < 0.25)) * cos(6*pi*z[,1])
beta2 <- function(z)  rep(1, nrow(z))
beta3 <- function(z)  10*sin(pi*z[,1]*z[,2]) + 20*(z[,3]-0.5)^2 + 10*z[,4] + 5*z[,5]

make_truth_list <- function(p) {
  stopifnot(p >= 3)
  truth <- vector("list", p + 1)
  truth[[1]] <- beta0
  truth[[2]] <- beta1
  truth[[3]] <- beta2
  truth[[4]] <- beta3
  if (p > 3) for (j in 4:p) truth[[j+1]] <- function(z) rep(0, nrow(z))
  truth
}


mse_beta_on_grid <- function(sum_obj, Zg, truth_list) {
  p1 <- dim(sum_obj)[3]; stopifnot(length(truth_list) == p1)
  mse_j <- numeric(p1)
  for (j in 1:p1) {
    truth_j <- truth_list[[j]](Zg[, 1:5, drop = FALSE])
    mse_j[j] <- mean((sum_obj[, "MEAN", j] - truth_j)^2)
  }
  list(per_j = mse_j, avg = mean(mse_j))
}
cov_beta_on_grid <- function(sum_obj, Zg, truth_list, level = 0.95) {
  p1 <- dim(sum_obj)[3]; stopifnot(length(truth_list) == p1)
  Lnm <- paste0("L", round(level*100)); Unm <- paste0("U", round(level*100))
  cov_j <- numeric(p1)
  for (j in 1:p1) {
    truth_j <- truth_list[[j]](Zg[, 1:5, drop = FALSE])
    L <- sum_obj[, Lnm, j]; U <- sum_obj[, Unm, j]
    cov_j[j] <- mean(truth_j >= L & truth_j <= U)
  }
  list(per_j = cov_j, avg = mean(cov_j))
}


## ---------- DGP (Sim 2) ----------
gen_data_sparse <- function(n_subj = 1000, m_per = 1, p = 50, R = 20, sigma = 1.0) {
  truth_list <- make_truth_list(p)
  N <- n_subj * m_per
  subj_id <- rep(seq_len(n_subj), each = m_per)
  Z <- matrix(runif(N * R), ncol = R)
  S <- outer(1:p, 1:p, function(a, b) 0.5^(abs(a - b)))
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = S)
  betas <- vapply(truth_list, function(f) f(Z[, 1:5, drop = FALSE]), numeric(N))
  xb <- rowSums(X * betas[, 2:(p + 1), drop = FALSE])
  Y <- betas[, 1] + xb + rnorm(N, 0, sigma)
  list(Y = Y, X = X, Z = Z, subj_id = subj_id, ni = rep(m_per, n_subj),
       truth_list = truth_list)
}

run_one_rep_sim2 <- function(seed,
                                        n_chains = 4,
                                        chain_seeds = NULL,
                                        n_subj = 1000, m_per = 1,
                                        p = 50, R = 20, sigma = 1.0,
                                        M_sparse = 50, M_van = 50,
                                        nd = 2000, burn = 400,
                                        n_cuts = 50, g = 200,
                                        slope_tau_mult = NULL,
                                        verbose = FALSE) {

  # Fix data & cutpoints once per replication
  set.seed(seed)
  dat <- gen_data_sparse(n_subj, m_per, p, R, sigma)
  cp  <- make_cutpoints(dat$Z, n_cuts = n_cuts)

  # Z-grid: vary z1,z2; others at 0.5
  Zg <- matrix(0.5, nrow = g, ncol = R)
  Zg[, 1] <- seq(0, 1, length.out = g)
  Zg[, 2] <- seq(0, 1, length.out = g)
  Zk_g <- matrix(0L, nrow = g, ncol = 1)

  if (is.null(chain_seeds)) {
    chain_seeds <- seed + 1000L * seq_len(n_chains)
  } else {
    stopifnot(length(chain_seeds) == n_chains)
  }

  tau_vec_sparse <- NULL
  if (!is.null(slope_tau_mult)) {
    tau0 <- 0.5 / sqrt(M_sparse)
    tau_vec_sparse <- c(tau0, rep(slope_tau_mult * tau0, p))
  }

  # run K chains
  per_chain <- do.call(rbind, lapply(seq_len(n_chains), function(k) {
    set.seed(chain_seeds[k])

    # --- sparseVCBART ---
    fit_sparse <- sparseVCBART::VCBART_ind(
      Y_train = dat$Y, subj_id_train = dat$subj_id, ni_train = dat$ni,
      X_train = dat$X,
      Z_cont_train = dat$Z, Z_cat_train = matrix(0L, 1, 1),
      cutpoints_list = cp,
      M = M_sparse, nd = nd, burn = burn,
      tau = tau_vec_sparse,
      save_samples = TRUE, save_trees = FALSE, verbose = verbose
    )

    # --- vanilla VCBART ---
    fit_van <- VCBART::VCBART_ind(
      Y_train = dat$Y, subj_id_train = dat$subj_id, ni_train = dat$ni,
      X_train = dat$X,
      Z_cont_train = dat$Z, Z_cat_train = matrix(0L, 1, 1),
      cutpoints_list = cp,
      tau = rep(0.5 / sqrt(M_van), p + 1),
      M = M_van, nd = nd, burn = burn,
      save_samples = TRUE, save_trees = FALSE, verbose = verbose
    )

    # predict betas on grid
    b_sparse_g <- sparseVCBART::predict_betas(fit_sparse, Z_cont = Zg, Z_cat = Zk_g, verbose = FALSE)
    b_van_g    <- VCBART::predict_betas(      fit_van,    Z_cont = Zg, Z_cat = Zk_g, verbose = FALSE)

    sum_sparse <- sparseVCBART::summarize_beta(b_sparse_g)
    sum_van    <- VCBART::summarize_beta(b_van_g)

    mse_s <- mse_beta_on_grid(sum_sparse, Zg, dat$truth_list)
    mse_v <- mse_beta_on_grid(sum_van,    Zg, dat$truth_list)
    cov_s <- cov_beta_on_grid(sum_sparse, Zg, dat$truth_list, level = 0.95)
    cov_v <- cov_beta_on_grid(sum_van,    Zg, dat$truth_list, level = 0.95)

    data.frame(
      chain = k,
      a_mse_sparse = mse_s$avg,
      a_mse_van    = mse_v$avg,
      b_cov_sparse = cov_s$avg,
      b_cov_van    = cov_v$avg,
      row.names = NULL
    )
  }))

  agg <- data.frame(
    seed = seed,
    a_mse_sparse = mean(per_chain$a_mse_sparse),
    a_mse_van    = mean(per_chain$a_mse_van),
    b_cov_sparse = mean(per_chain$b_cov_sparse),
    b_cov_van    = mean(per_chain$b_cov_van)
  )

  agg$a_mse_sparse_se <- stats::sd(per_chain$a_mse_sparse)/sqrt(n_chains)
  agg$a_mse_van_se    <- stats::sd(per_chain$a_mse_van)/sqrt(n_chains)
  agg$b_cov_sparse_se <- stats::sd(per_chain$b_cov_sparse)/sqrt(n_chains)
  agg$b_cov_van_se    <- stats::sd(per_chain$b_cov_van)/sqrt(n_chains)

  list(aggregate = agg, per_chain = per_chain)
}

out <- run_one_rep_sim2(
  seed = 1, n_chains = 4,
  n_subj = 1000, m_per = 1, p = 50, R = 20, sigma = 1.0,
  M_sparse = 50, M_van = 50, nd = 2000, burn = 400,
  n_cuts = 50, g = 200, slope_tau_mult = NULL, verbose = FALSE
)

#out$aggregate   
