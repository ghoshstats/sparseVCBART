## ===================== Simulation 2 (VCBART variants) ======================
## Intercept + first 3 slopes as in Sim-1; all remaining slopes are zero.
suppressPackageStartupMessages({
  library(MASS)
  library(VCBART)
  library(sparseVCBART)
})
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
stopifnot(!is.na(seed), seed >= 0)

rmse_fun <- function(sum_obj, j, truth_fun, Zg) {
  mean((sum_obj[, "MEAN", j] - truth_fun(Zg[, 1:2, drop = FALSE]))^2)^0.5
}
width95 <- function(sum_obj, j) mean(sum_obj[, "U95", j] - sum_obj[, "L95", j])


## ---------- Helpers ----------
make_cutpoints <- function(Z_cont, n_cuts = 50) {
  if (!is.matrix(Z_cont)) Z_cont <- as.matrix(Z_cont)
  out <- vector("list", ncol(Z_cont))
  for (r in seq_len(ncol(Z_cont))) {
    z  <- Z_cont[, r]
    qs <- unique(as.numeric(stats::quantile(z, probs = seq(0, 1, length.out = n_cuts + 2))))
    out[[r]] <- qs[-c(1, length(qs))]
  }
  out
}

# Sim-1 truths
beta0 <- function(z)  3*z[,1] + ((2 - 5*(z[,2] > 0.5)) * sin(pi*z[,1])) - 2*(z[,2] > 0.5)
beta1 <- function(z)  ((3 - 3*z[,2]) * (z[,1] > 0.6) - 10*sqrt(pmax(z[,1], 0)) * (z[,1] < 0.25)) * cos(6*pi*z[,1])
beta2 <- function(z)  rep(1, nrow(z))
beta3 <- function(z)  10*sin(pi*z[,1]*z[,2]) + 20*(z[,3]-0.5)^2 + 10*z[,4] + 5*z[,5]

# Extend to p slopes: keep beta0..beta3, set the rest to zero functions
make_truth_list <- function(p) {
  stopifnot(p >= 3)
  truth <- vector("list", p + 1)
  truth[[1]] <- beta0
  truth[[2]] <- beta1
  truth[[3]] <- beta2
  truth[[4]] <- beta3
  if (p > 3) {
    for (j in 4:p) truth[[j+1]] <- function(z) rep(0, nrow(z))
  }
  truth
}

# Metrics on a Z-grid
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

## ---------- DGP ----------
gen_data_sparse_like_sim1 <- function(n_subj = 1000, m_per = 1, p = 50, R = 20, sigma = 1.0) {
  truth_list <- make_truth_list(p)
  N <- n_subj * m_per
  subj_id <- rep(seq_len(n_subj), each = m_per)
  Z <- matrix(runif(N * R), ncol = R)
  # same correlation as Sim-1
  S <- outer(1:p, 1:p, function(a, b) 0.5^(abs(a - b)))
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = S)
  betas <- vapply(truth_list, function(f) f(Z[, 1:5, drop = FALSE]), numeric(N))
  xb <- rowSums(X * betas[, 2:(p + 1), drop = FALSE])
  Y <- betas[, 1] + xb + rnorm(N, 0, sigma)
  list(Y = Y, X = X, Z = Z, subj_id = subj_id, ni = rep(m_per, n_subj),
       truth_list = truth_list)
}

## ---------- One replication (a,b only) ----------
run_one_rep_sim2_ab <- function(seed,
                                n_subj = 1000, m_per = 1,
                                p = 50, R = 20, sigma = 1.0,
                                M_sparse = 50, M_van = 50,
                                nd = 2000, burn = 400,
                                n_cuts = 50, g = 200,
                                slope_tau_mult = NULL,
                                verbose = FALSE) {
  set.seed(seed)
  dat <- gen_data_sparse_like_sim1(n_subj, m_per, p, R, sigma)
  
  # cutpoints from full Z
  cp <- make_cutpoints(dat$Z, n_cuts = n_cuts)
  
  # sparseVCBART: optionally inflate slope taus
  tau_vec_sparse <- NULL
  if (!is.null(slope_tau_mult)) {
    tau0 <- 0.5 / sqrt(M_sparse)
    tau_vec_sparse <- c(tau0, rep(slope_tau_mult * tau0, p))
  }
  
  fit_sparse <- sparseVCBART::VCBART_ind(
    Y_train = dat$Y,
    subj_id_train = dat$subj_id,
    ni_train = dat$ni,
    X_train = dat$X,
    Z_cont_train = dat$Z,
    Z_cat_train = matrix(0L, 1, 1),
    cutpoints_list = cp,
    M = M_sparse, nd = nd, burn = burn,
    tau = tau_vec_sparse,
    save_samples = TRUE, save_trees = TRUE, verbose = verbose
  )
  
  fit_van <- VCBART::VCBART_ind(
    Y_train = dat$Y,
    subj_id_train = dat$subj_id,
    ni_train = dat$ni,
    X_train = dat$X,
    Z_cont_train = dat$Z,
    Z_cat_train = matrix(0L, 1, 1),
    cutpoints_list = cp,
    tau = rep(0.5 / sqrt(M_van), p + 1),
    M = M_van, nd = nd, burn = burn,
    save_samples = TRUE, save_trees = TRUE, verbose = verbose
  )
  
  ## Z-grid: vary Z1 and Z2, hold others at 0.5 (as in Sim-1)
  Zg <- matrix(0.5, nrow = g, ncol = R)
  Zg[, 1] <- seq(0, 1, length.out = g)
  Zg[, 2] <- seq(0, 1, length.out = g)
  
  b_sparse_g <- sparseVCBART::predict_betas(fit_sparse, Z_cont = Zg, Z_cat = matrix(0L, g, 1))
  b_van_g    <- VCBART::predict_betas(      fit_van,    Z_cont = Zg, Z_cat = matrix(0L, g, 1))
  
  sum_sparse <- sparseVCBART::summarize_beta(b_sparse_g)   # [g x {MEAN,L95,U95} x (p+1)]
  sum_van    <- VCBART::summarize_beta(b_van_g)
  
  # metrics (average across j = 0..p and grid)
  mse_s <- mse_beta_on_grid(sum_sparse, Zg, dat$truth_list)
  mse_v <- mse_beta_on_grid(sum_van,    Zg, dat$truth_list)
  cov_s <- cov_beta_on_grid(sum_sparse, Zg, dat$truth_list, level = 0.95)
  cov_v <- cov_beta_on_grid(sum_van,    Zg, dat$truth_list, level = 0.95)
  
  data.frame(
    seed = seed,
    a_mse_sparse = mse_s$avg,
    a_mse_van    = mse_v$avg,
    b_cov_sparse = cov_s$avg,
    b_cov_van    = cov_v$avg
  )
}

res <- run_one_rep_sim2_ab(seed)
print(res)
saveRDS(res, file = sprintf("new_sim2_sparssevanillavcbart_seed%02d.rds", seed))

## ===================== Example run ===================== 
# (Adjust p, nd/burn, etc. as needed)
# set.seed(42)
# res_ab <- run_many_sim2_ab(
#   seeds = 1:25,
#   n_subj = 1000, m_per = 1,
#   p = 50, R = 20, sigma = 1.0,
#   M_sparse = 50, M_van = 50,
#   nd = 2000, burn = 400,
#   n_cuts = 50, g = 200,
#   slope_tau_mult = NULL,   # or e.g., 1.5 or 2.0 if you want looser sparse slopes
#   verbose = FALSE
# )
# print(res_ab)
# colMeans(res_ab[ , -1, drop = FALSE])  # quick averages across seeds
# saveRDS(res_ab, file = "sim2_variant_ab_results.rds")
