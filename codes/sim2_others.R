## ======================== Sim-2 settings for other methods ========================
## Truth: beta0..beta3 as in Sim-1; all remaining slopes are zero functions.
## Metrics: (a) average MSE for beta_j(z); (b) average 95% coverage for beta_j(z)
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1]); stopifnot(!is.na(seed), seed >= 0)

library(MASS)
source("competitor_wrappers.R")  

## ---------- Sim-1 truths ----------
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

## ---------- DGP (Sim 2) ----------
gen_data_sim2 <- function(n_subj=1000, m_per=1, p=50, R=20, sigma=1.0) {
  N <- n_subj * m_per
  Z <- matrix(runif(N * R), ncol = R)
  S <- outer(1:p, 1:p, function(a,b) 0.5^(abs(a-b)))   # SAME as Sim-1
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = S)
  truth_list <- make_truth_list(p)
  betas <- vapply(truth_list, function(f) f(Z[, 1:5, drop = FALSE]), numeric(N))
  Y <- betas[,1] + rowSums(X * betas[, 2:(p+1), drop=FALSE]) + rnorm(N, 0, sigma)
  list(Y=Y, X=X, Z=Z, truth_list=truth_list)
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

#################### Run BTVCM ########################
run_one_rep_sim2_btvcm <- function(seed,
                                      n_subj = 1000, p = 50, R = 20, sigma = 1.0,
                                      g = 200,
                                      lambda = 0.05, ntree = 250, B = 40) {
  set.seed(seed)
  dat <- gen_data_sim2(n_subj=n_subj, m_per=1, p=p, R=R, sigma=sigma)
  
  ## Z-grid: vary Z1 & Z2; others fixed at 0.5 (truths use first 5 modifiers)
  Zg <- matrix(0.5, nrow=g, ncol=R)
  Zg[,1] <- seq(0, 1, length.out=g)
  Zg[,2] <- seq(0, 1, length.out=g)
  
  ## Ask boosted TVCM for beta surfaces on the grid: pass X_test = 0
  Xg0 <- matrix(0, nrow=g, ncol=p)
  res <- boosted_tvcm_wrapper(
    Y_train = dat$Y,
    X_train = dat$X,
    Z_train = dat$Z,
    X_test  = Xg0,
    Z_test  = Zg,
    intercept = TRUE,
    lambda = lambda,
    ntree = ntree,
    B = B           # must be > 0 to get intervals
  )
  if (is.null(res) || is.null(res$test$beta))
    stop(sprintf("btvcm wrapper failed or returned no beta for seed=%d.", seed))
  
  bst <- res$test$beta  # expected dims: [g x {MEAN,SD,L95,U95} x (p+1)]
  need <- c("MEAN","L95","U95")
  have <- dimnames(bst)[[2]]
  stopifnot(all(need %in% have))
  # Conform to [g x {MEAN,L95,U95} x (p+1)]
  ss <- array(NA_real_, dim = c(dim(bst)[1], 3, dim(bst)[3]),
              dimnames = list(NULL, c("MEAN","L95","U95"), NULL))
  ss[, "MEAN", ] <- bst[, "MEAN", ]
  ss[, "L95",  ] <- bst[, "L95",  ]
  ss[, "U95",  ] <- bst[, "U95",  ]
  
  # Metrics (average across j=0..p and grid)
  mse <- mse_beta_on_grid(ss, Zg, dat$truth_list)
  cov <- cov_beta_on_grid(ss, Zg, dat$truth_list, level = 0.95)
  
  data.frame(
    method = "BTVMC",
    seed   = seed,
    a_mse_beta = mse$avg,
    b_cov_beta = cov$avg,
    stringsAsFactors = FALSE
  )
}
res_BTVCM <- run_one_rep_sim2_btvcm(seed,B=50)
#print(res_BTVCM)

#################### Run LM ########################
run_one_rep_sim2_LM <- function(seed,
                                        n_subj = 1000, m_per = 1,
                                        p = 50, R = 20, sigma = 1,
                                        d0 = 5, s0 = 2,
                                        g = 200) {
  set.seed(seed)
  
  tr <- gen_data_hd(n_subj = n_subj, m_per = m_per, p = p, R = R,
                    sigma = sigma, d0 = d0, s0 = s0)
  
  Zg <- matrix(0.5, nrow = g, ncol = R); Zg[, 1] <- seq(0, 1, length.out = g)
  
  Xg0 <- matrix(0, nrow = g, ncol = p)
  
  fit <- lm_wrapper(tr$Y, tr$X, tr$Z, Xg0, Zg)
  
  ss2 <- fit$test$beta                    # dims: g x 4 x (p+1)
  stopifnot(all(c("MEAN","L95","U95") %in% dimnames(ss2)[[2]]))
  ss2 <- ss2[, c("MEAN","L95","U95"), , drop = FALSE]  # ensure right order
  
  # identify active vs zero (intercept + first d0 are active)
  is_zero <- c(FALSE, rep(FALSE, d0), rep(TRUE, p - d0))  # length p+1
  
  truth_list <- gen_sparse_truth(p = p, R = R, d0 = d0, s0 = s0)
  
  rmse_act <- rmse_zero <- wid_act <- wid_zero <- numeric(0)
  for (j in 1:(p + 1)) {
    f <- truth_list[[j]]
    if (is_zero[j]) {
      rmse_zero <- c(rmse_zero, rmse_fun(ss2, j, f, Zg))
      wid_zero  <- c(wid_zero,  width95(ss2, j))
    } else {
      rmse_act <- c(rmse_act, rmse_fun(ss2, j, f, Zg))
      wid_act  <- c(wid_act,  width95(ss2, j))
    }
  }
  
  data.frame(
    method = "LM",
    seed   = seed,
    rmse_act_mean   = mean(rmse_act),
    rmse_zero_mean  = mean(rmse_zero),
    width_act_mean  = mean(wid_act),
    width_zero_mean = mean(wid_zero),
    row.names = NULL
  )
}

# lm_actzero_25 <- do.call(rbind, lapply(1:25, function(s)
#   run_one_rep_sim2_LM(seed = s))
# )
