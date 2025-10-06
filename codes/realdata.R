######################### sparseVCBART real data analysis ########################
##################################################################################
load("e1_data_for_wrapper.RData")
library(dplyr)
library(ggplot2)
library(tidyr)
library(rpart)
library(rpart.plot)

rmse <- function(yhat, y) sqrt(mean((yhat - y)^2))
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
make_cutpoints <- function(Z_cont, n_cuts = 50) {
  if (!is.matrix(Z_cont)) Z_cont <- as.matrix(Z_cont)
  out <- vector("list", ncol(Z_cont))
  for (r in seq_len(ncol(Z_cont))) {
    z <- Z_cont[, r]
    qs <- unique(as.numeric(stats::quantile(z, probs = seq(0, 1, length.out = n_cuts + 2))))
    out[[r]] <- qs[-c(1, length(qs))]
  }
  out
}

make_cat_levels <- function(Z_cat) {
  if (!is.matrix(Z_cat)) Z_cat <- as.matrix(Z_cat)
  lapply(seq_len(ncol(Z_cat)), function(j) sort(unique(Z_cat[, j])))
}

## Combine posterior mean betas into posterior mean of E[Y | X, Z]
## sum_obj: output from summarize_beta(...): [n_test x {MEAN,L95,U95} x (p+1)]
yhat_from_beta_summ <- function(sum_obj, X_test) {
  beta0 <- sum_obj[, "MEAN", 1]
  if (ncol(X_test) == 0) return(beta0)
  B     <- sapply(1:ncol(X_test), function(j) sum_obj[, "MEAN", j + 1])  # [n_test x p]
  as.numeric(beta0 + rowSums(X_test * B))
}

######################## Adding noise variables #####################

## --------------------------
## Add 18 noise predictors to X (p = 20 total)
## --------------------------
set.seed(123)
stopifnot(nrow(X1_all) == length(Y_all),
          nrow(Z_cat_all) == length(Y_all),
          nrow(Z_cont_all) == length(Y_all))

n  <- nrow(X1_all)
p0 <- ncol(X1_all)       
p  <- 20                 

# Create 18 independent noise columns (standard normal)
p_noise <- p - p0
X_noise <- matrix(rnorm(n * p_noise), nrow = n, ncol = p_noise)
colnames(X_noise) <- paste0("noise_", seq_len(p_noise))

# Augmented design: first 2 are the true predictors
X_aug <- cbind(X1_all, X_noise)
colnames(X_aug)[1:p0] <- colnames(X1_all)

# Build if missing
if (!exists("cp_list")) cp_list <- make_cutpoints(Z_cont_all, n_cuts = 50)
if (!exists("cat_levels_list")) cat_levels_list <- make_cat_levels(Z_cat_all)

## --------------------------
## Fit sparseVCBART on augmented X
## --------------------------
M    <- 50
nd   <- 1200
burn <- 400

# VCBART default tau_j per paper: 0.5/sqrt(M)
tau_vec <- rep(0.5 / sqrt(M), ncol(X_aug) + 1)  # +1 for intercept

fit_sparse <- sparseVCBART::VCBART_ind(
  Y_train        = Y_all,
  subj_id_train  = subj_id_all,
  ni_train       = ni_all,
  X_train        = X_aug,
  Z_cont_train   = Z_cont_all,
  Z_cat_train    = Z_cat_all,
  cutpoints_list = cp_list,
  cat_levels_list = cat_levels_list,
  M = M, nd = nd, burn = burn, 
  save_samples = TRUE, save_trees = TRUE, verbose = FALSE
)

## --------------------------
## Build a realistic Z grid to evaluate beta_j(z)
## --------------------------
set.seed(456)
g <- min(300L, n)  # grid size
idx_grid     <- sample.int(n, g)
Z_cont_grid  <- as.matrix(Z_cont_all[idx_grid, , drop = FALSE])
Z_cat_grid   <- as.matrix(Z_cat_all[idx_grid, , drop = FALSE])

# Predict betas on grid
b_sparse_draws <- sparseVCBART::predict_betas(
  fit_sparse, Z_cont = Z_cont_grid, Z_cat = Z_cat_grid, verbose = FALSE
)
sum_sparse <- sparseVCBART::summarize_beta(b_sparse_draws)

## --------------------------
## Selection scores and top-2 selection (exclude intercept)
## --------------------------
# sparseVCBART: score = median(lambda_j)
lam_draws <- fit_sparse$lambda        # [iter x (p+1)], includes intercept in column 1
stopifnot(ncol(lam_draws) == ncol(X_aug) + 1)
lam_med   <- apply(lam_draws, 2, median)
lam_med_pred <- lam_med[-1]         
names(lam_med_pred) <- colnames(X_aug)

true_idx  <- 1:p0               
zero_idx  <- (p0+1):p
sel_sparse_id <- match(sel_sparse_vars, colnames(X_aug))
print(sel_sparse_id)
                  
####################### Fit-the-fit strategy ###################################

## Posterior mean beta_j(z) at observed modifiers -------------------------

# Predict beta surfaces at the observed Z (no need for X)
b_draws <- sparseVCBART::predict_betas(
  fit_sparse,
  Z_cont = Z_cont_all,
  Z_cat  = Z_cat_all,
  verbose = FALSE
)
b_summ  <- sparseVCBART::summarize_beta(b_draws)

beta_mean_df <- as.data.frame(b_summ[, "MEAN", ])
colnames(beta_mean_df) <- c("(Intercept)", colnames(X1_all))

##  Which modifiers did sparseVCBART actually use? -------------------------

vc_summarize <- function(varcounts) {
  dm <- dim(varcounts)
  if (length(dm) == 3L) {
    # assume [iter x R x (p+1)]
    iter <- dm[1]; Rtot <- dm[2]; p1 <- dm[3]
    vc_mean <- apply(varcounts, c(2,3), mean, na.rm = TRUE)
    vc_prob <- apply(varcounts > 0, c(2,3), mean, na.rm = TRUE)
  } else if (length(dm) == 2L) {
    # single matrix: treat as deterministic counts; prob = 1 if >0 else 0
    Rtot <- dm[1]; p1 <- dm[2]
    vc_mean <- varcounts
    vc_prob <- (varcounts > 0) * 1
  } else {
    stop("Unrecognized shape for fit_sparse$varcounts")
  }
  list(mean = vc_mean, prob = vc_prob)
}

vc <- vc_summarize(fit_sparse$varcounts)

Z_cont_df <- as.data.frame(Z_cont_all)
Z_cat_df  <- as.data.frame(Z_cat_all)
Z_cat_df[] <- lapply(Z_cat_df, function(x) factor(x, levels = sort(unique(x))))
Z_all_df  <- cbind(Z_cont_df, Z_cat_df)
mod_names <- colnames(Z_all_df)
stopifnot(nrow(Z_all_df) == nrow(beta_mean_df))


# Use posterior probability-of-use > 0.25 by default.
select_modifiers <- function(j_col, prob_mat, mean_mat, names_all,
                             prob_thresh = 0.25, top_k_fallback = 3) {
  pr <- prob_mat[, j_col]
  keep <- which(pr > prob_thresh)
  if (length(keep) == 0L) {
    # fallback: top-k by expected split count
    ord <- order(mean_mat[, j_col], decreasing = TRUE, na.last = NA)
    keep <- head(ord, top_k_fallback)
  }
  names_all[keep]
}

p <- ncol(X1_all)                              
p1 <- p + 1                                   
sel_list <- vector("list", length = p1)
names(sel_list) <- colnames(beta_mean_df)

for (j in 1:p1) {
  sel_list[[j]] <- select_modifiers(j, vc$prob, vc$mean, mod_names,
                                    prob_thresh = 0.25, top_k_fallback = 3)
}

## Fit CARTs to posterior means using only selected modifiers --------------

fit_tree_for_beta <- function(y, Z_df, sel_names,
                              minsplit = 40, cp = 0.0, xval = 10) {
  if (length(sel_names) == 0L) return(NULL)
  dat <- data.frame(y = y, Z_df[, sel_names, drop = FALSE])
  fit <- rpart(y ~ ., data = dat, method = "anova",
               control = rpart.control(minsplit = minsplit, cp = cp, xval = xval))
  cpt <- printcp(fit)
  best_row <- which.min(cpt[,"xerror"])
  se_rule <- cpt[best_row, "xerror"] + cpt[best_row, "xstd"]
  cp_1se  <- max(cpt[cpt[,"xerror"] <= se_rule, "CP"])
  pruned  <- prune(fit, cp = cp_1se)
  list(full = fit, pruned = pruned, cp_table = cpt)
}

trees <- vector("list", length = p1)
names(trees) <- colnames(beta_mean_df)

for (j in 1:p1) {
  yj <- beta_mean_df[[j]]           # posterior mean beta_j(Z)
  trees[[j]] <- fit_tree_for_beta(y = yj, Z_df = Z_all_df, sel_names = sel_list[[j]])
}

##  How good is the “fit-the-fit”?  ----------------------------------------

tree_rmse <- function(tree_obj, y, Z_df, sel_names) {
  if (is.null(tree_obj)) return(NA_real_)
  dat <- data.frame(Z_df[, sel_names, drop = FALSE])
  pred <- predict(tree_obj, newdata = dat)
  sqrt(mean((y - pred)^2, na.rm = TRUE))
}

pruned_list <- lapply(trees, `[[`, "pruned")          
y_list      <- as.list(beta_mean_df)                  
sel_per_j   <- sel_list                              

stopifnot(identical(names(pruned_list), names(y_list)),
          identical(names(pruned_list), names(sel_per_j)))

rmse_vals <- mapply(
  FUN      = tree_rmse,
  tree_obj = pruned_list,
  y        = y_list,
  sel_names= sel_per_j,
  MoreArgs = list(Z_df = Z_all_df),
  SIMPLIFY = TRUE
)

rmse_table <- data.frame(
  coefficient = names(trees),
  n_modifiers = vapply(sel_list, length, integer(1)),
  rmse_pruned = as.numeric(rmse_vals),
  stringsAsFactors = FALSE
)

print(rmse_table)
# if (!is.null(trees[["treat_pg_apt"]])) {
#   rpart.plot(trees[["treat_pg_apt"]][["pruned"]],
#              type = 2, extra = 101, under = TRUE, faclen = 0,
#              main = "")
# }


