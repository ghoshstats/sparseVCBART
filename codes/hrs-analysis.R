
library(sparseVCBART)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)


M       <- 50
nd      <- 2000
burn    <- 400
thin    <- 1
n_chains <- 4
base_seed <- 424

load("hrs_data1.RData")

subj_id_all <- rep(seq_len(n_all), times = ni_all)

stopifnot(
  length(Y_all) == nrow(X1_all),
  length(Y_all) == nrow(Z_cont_all),
  length(Y_all) == nrow(Z_cat_all),
  length(Y_all) == length(subj_id_all)
)

keep_idx <- function(n_iter, burn) {
  stopifnot(burn < n_iter)
  (burn + 1):n_iter
}

drop_burn_array <- function(arr, burn) {
  dm <- dim(arr)
  if (is.null(dm)) stop("Object has no dim()")
  n_iter <- dm[1]
  idx <- keep_idx(n_iter, burn)
  if (length(dm) == 2) return(arr[idx, , drop = FALSE])
  if (length(dm) == 3) return(arr[idx, , , drop = FALSE])
  if (length(dm) == 4) return(arr[idx, , , , drop = FALSE])
  stop("Unsupported dimension for array")
}

run_one_chain <- function(chain_num, X_train, tag) {
  set.seed(base_seed + chain_num)
  p <- ncol(X_train)
  mu0 <- rep(0, p + 1)
  tau <- rep(0.5 / sqrt(M), p + 1)
  
  fit <- sparseVCBART::VCBART_cs(
    Y_train         = Y_all,
    subj_id_train   = subj_id_all,
    ni_train        = ni_all,
    X_train         = X_train,
    Z_cont_train    = Z_cont_all,
    Z_cat_train     = Z_cat_all,
    unif_cuts       = unif_cuts,
    cutpoints_list  = cutpoints_list,
    cat_levels_list = cat_levels_list,
    sparse = TRUE,
    M = M, mu0 = mu0, tau = tau,
    nd = nd, burn = burn, thin = thin,
    save_trees = TRUE, save_samples = TRUE,
    verbose = TRUE
  )
  
  lam_post <- {
    n_iter <- nrow(fit$lambda)
    fit$lambda[keep_idx(n_iter, burn), , drop = FALSE]
  }
  
  vc_post <- drop_burn_array(fit$varcounts, burn)
  
  chain_obj <- list(
    chain_num = chain_num,
    tag = tag,
    lambda_post = lam_post,
    varcounts_post = vc_post,
    sigma = fit$sigma
  )
  
  saveRDS(fit, file = sprintf("HRS_%s_chain%d_fit.rds", tag, chain_num))
  save(chain_obj, file = sprintf("HRS_%s_chain%d.RData", tag, chain_num))
  
  invisible(chain_obj)
}

chains_full <- vector("list", n_chains)
for (cc in 1:n_chains) {
  chains_full[[cc]] <- run_one_chain(chain_num = cc, X_train = X1_all, tag = "full")
}

lam_all_full <- do.call(rbind, lapply(chains_full, `[[`, "lambda_post"))
lam_med_full <- apply(lam_all_full, 2, median)
lam_med_full_pred <- lam_med_full[-1]  # drop intercept
names(lam_med_full_pred) <- colnames(X1_all)
lam_med_full_pred

drop_vars <- c("INCOME", "LBR")
stopifnot(all(drop_vars %in% colnames(X1_all)))

X1_reduced <- X1_all[, !(colnames(X1_all) %in% drop_vars), drop = FALSE]

chains_red <- vector("list", n_chains)
for (cc in 1:n_chains) {
  chains_red[[cc]] <- run_one_chain(chain_num = cc, X_train = X1_reduced, tag = "reduced")
}

lam_all_red <- do.call(rbind, lapply(chains_red, `[[`, "lambda_post"))
lam_med_red <- apply(lam_all_red, 2, median)
lam_med_red_pred <- lam_med_red[-1]
names(lam_med_red_pred) <- colnames(X1_reduced)
lam_med_red_pred


############################# Lambda plots for full and reduced models ################################
#######################################################################################################
fit_full <- readRDS("HRS_sparseVCBART.rds")
fit_red  <- readRDS("HRS_sparseVCBART_noINCOME_noLBR.rds")

tidy_lambda <- function(fit, pred_names, scenario, drop_intercept = TRUE) {
  lam <- fit$lambda
  if (drop_intercept) lam <- lam[, -1, drop = FALSE]
  colnames(lam) <- pred_names
  
  as.data.frame(lam) %>%
    pivot_longer(everything(), names_to = "covariate", values_to = "lambda") %>%
    mutate(scenario = scenario)
}

df <- bind_rows(df_full, df_red)

cov_order <- c("CHLD_SES","HS_COMPL","COLLEGE_COMPL","INCOME","LBR")
df <- df %>%
  mutate(covariate = factor(covariate, levels = cov_order))

tidy_lambda <- function(fit, pred_names, scenario) {
  lam <- fit$lambda[, -1, drop = FALSE]   # drop intercept
  colnames(lam) <- pred_names
  
  as.data.frame(lam) %>%
    pivot_longer(everything(), names_to = "covariate", values_to = "lambda") %>%
    mutate(scenario = scenario)
}

df_full <- tidy_lambda(fit_full, colnames(X1_all),     "Expanded Model (p = 5)")
df_red  <- tidy_lambda(fit_red,  colnames(X1_reduced), "Initial Model (p = 3)")
df <- bind_rows(df_full, df_red)

df <- df %>%
  mutate(covariate = factor(covariate, levels = cov_order))

sig_thresh <- 10

sig_flag <- df %>%
  group_by(scenario, covariate) %>%
  summarize(lam_med = median(lambda, na.rm = TRUE), .groups = "drop") %>%
  mutate(sig = lam_med >= sig_thresh)

df_plot <- df %>%
  left_join(sig_flag %>% select(scenario, covariate, sig),
            by = c("scenario", "covariate")) %>%
  mutate(sig_lab = if_else(sig, "significant", "insignificant"))

df_plot <- df_plot %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Initial Model (p = 3)", "Expanded Model (p = 5)")
    )
  )

p_box <- ggplot(df_plot, aes(x = covariate, y = lambda, fill = sig_lab)) +
  geom_boxplot(width = 0.65, outlier.alpha = 0.10) +
  geom_hline(yintercept = sig_thresh, linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = c("significant" = "red", "insignificant" = "grey70")) +
  labs(
    x = NULL,
    y = expression(lambda[j] * " (log10 scale)"),
    fill = NULL,
    title = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "bottom",
    
    panel.border = element_rect(fill = NA, linewidth = 0.9),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(1.0, "lines"),
    
    strip.background = element_rect(fill = "grey95", linewidth = 0.6),
    strip.placement = "outside"
  )

p_box
ggsave("lambda_full_vs_reduced.png", p_box, width = 7, height = 4.2, dpi = 300)

