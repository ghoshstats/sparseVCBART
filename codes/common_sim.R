## ===== Common helpers for all methods =====
library(MASS)

## ---------- Helper: cutpoints ----------
make_cutpoints <- function(Z_cont, n_cuts = 25) {
  if (!is.matrix(Z_cont)) Z_cont <- as.matrix(Z_cont)
  out <- vector("list", ncol(Z_cont))
  for (r in seq_len(ncol(Z_cont))) {
    z  <- Z_cont[, r]
    qs <- unique(as.numeric(stats::quantile(z, probs = seq(0, 1, length.out = n_cuts + 2))))
    out[[r]] <- qs[-c(1, length(qs))]
  }
  out
}

## ---------- Truth -----------------
beta0 <- function(z)  3*z[,1] + ( (2 - 5*(z[,2] > 0.5)) * sin(pi*z[,1]) ) - 2*(z[,2] > 0.5)
beta1 <- function(z)  ((3 - 3*z[,2]) * (z[,1] > 0.6) - 10*sqrt(pmax(z[,1], 0)) * (z[,1] < 0.25)) * cos(6*pi*z[,1])
beta2 <- function(z)  rep(1, nrow(z))
beta3 <- function(z)  10*sin(pi*z[,1]*z[,2]) + 20*(z[,3]-0.5)^2 + 10*z[,4] + 5*z[,5]
truth_list <- list(beta0, beta1, beta2, beta3)

## ------------- DGP ------------------
gen_data <- function(n_subj = 1000, m_per = 1, p = 3, R = 20, sigma = 1.0) {
  N <- n_subj * m_per
  subj_id <- rep(seq_len(n_subj), each = m_per)
  Z <- matrix(runif(N * R), ncol = R)
  S <- outer(1:p, 1:p, function(a, b) 0.5^(abs(a - b)))
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = S)
  betas <- vapply(truth_list, function(f) f(Z[, 1:5, drop = FALSE]), numeric(N))
  xb <- rowSums(X * betas[, 2:(p + 1), drop = FALSE])
  Y <- betas[, 1] + xb + rnorm(N, 0, sigma)
  list(Y = Y, X = X, Z = Z, subj_id = subj_id, ni = rep(m_per, n_subj), betas = betas)
}

## ---------- Metrics ------------------
# Mean prediction from beta summaries
yhat_from_beta_summ <- function(summ, X) {
  n  <- dim(summ)[1]; p1 <- dim(summ)[3]; p <- p1 - 1
  mu <- summ[, "MEAN", 1]
  if (p > 0) for (j in 1:p) mu <- mu + summ[, "MEAN", j + 1] * X[, j]
  mu
}

# MSE over a grid Zg against truth function list
mse_beta_on_grid <- function(sum_obj, Zg, truth_list) {
  p1 <- dim(sum_obj)[3]; stopifnot(length(truth_list) == p1)
  mse_j <- numeric(p1)
  for (j in 1:p1) {
    truth_j <- truth_list[[j]](Zg[, 1:5, drop = FALSE])
    mse_j[j] <- mean((sum_obj[, "MEAN", j] - truth_j)^2)
  }
  list(per_j = mse_j, avg = mean(mse_j))
}

# Coverage for beta intervals on grid
cov_beta_on_grid <- function(sum_obj, Zg, truth_list, level = 0.95) {
  p1 <- dim(sum_obj)[3]; stopifnot(length(truth_list) == p1)
  Lname <- paste0("L", round(level*100)); Uname <- paste0("U", round(level*100))
  cov_j <- numeric(p1)
  for (j in 1:p1) {
    truth_j <- truth_list[[j]](Zg[, 1:5, drop = FALSE])
    L <- sum_obj[, Lname, j]; U <- sum_obj[, Uname, j]
    cov_j[j] <- mean(truth_j >= L & truth_j <= U)
  }
  list(per_j = cov_j, avg = mean(cov_j))
}

# Predictive metrics via posterior predictive draws (for VCBART-variants)
pred_metrics <- function(beta_draws, sigma_draws, X_te, Y_te, level = 0.95) {
  nd <- dim(beta_draws)[1]; n <- dim(beta_draws)[2]; p1 <- dim(beta_draws)[3]; p <- p1 - 1
  sigma_draws <- as.numeric(sigma_draws[seq_len(nd)])
  mu_mean <- rep(0, n)
  z_draws <- matrix(rnorm(nd * n, 0, 1), nd, n)
  yrep <- matrix(NA_real_, nd, n)
  for (d in 1:nd) {
    mu_d <- beta_draws[d, , 1]
    if (p > 0) for (j in 1:p) mu_d <- mu_d + beta_draws[d, , j + 1] * X_te[, j]
    mu_mean <- mu_mean + mu_d / nd
    yrep[d, ] <- mu_d + sigma_draws[d] * z_draws[d, ]
  }
  alpha <- (1 - level) / 2
  L <- apply(yrep, 2, stats::quantile, probs = alpha, names = FALSE)
  U <- apply(yrep, 2, stats::quantile, probs = 1 - alpha, names = FALSE)
  cov_pred <- mean(Y_te >= L & Y_te <= U)
  rmse <- sqrt(mean((mu_mean - Y_te)^2))
  list(rmse = rmse, cov = cov_pred)
}

pred_from_wrapper <- function(ystar_mat, Y_te) {
  if (is.null(ystar_mat)) return(c(rmse = NA_real_, cov = NA_real_))
  mu <- ystar_mat[, "MEAN"]; L <- ystar_mat[, "L95"]; U <- ystar_mat[, "U95"]
  c(rmse = sqrt(mean((mu - Y_te)^2)), cov = mean(Y_te >= L & Y_te <= U))
}
