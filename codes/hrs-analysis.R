#args <- commandArgs(TRUE)
#chain_num <- as.numeric(args[1])
library(sparseVCBART)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

chain_num <- 1
set.seed(424+chain_num)
load("hrs_data1.RData")

p <- ncol(X1_all)
R_cont <- ncol(Z_cont_all)

R_cat <- ncol(Z_cat_all)

R <- R_cont + R_cat

M <- 50
mu0 <- rep(0, times = p+1)
tau <- rep(0.5/sqrt(M), times = p+1)

subj_id_all <- c()
for(i in 1:n_all){
  subj_id_all <- c(subj_id_all, rep(i, times = ni_all[i]))
}

nd <- 1000
burn <- 1000
thin <- 1

fit <- sparseVCBART::VCBART_cs(Y_train = Y_all,
                         subj_id_train = subj_id_all,
                         ni_train = ni_all,
                         X_train = X1_all,
                         Z_cont_train = Z_cont_all,
                         Z_cat_train = Z_cat_all,
                         unif_cuts = unif_cuts,
                         cutpoints_list = cutpoints_list,
                         cat_levels_list = cat_levels_list,
                         sparse = TRUE,
                         M = 50, mu0 = mu0, tau = tau,
                         nd = nd, burn = burn, thin = thin,
                         save_trees = TRUE, save_samples = TRUE,
                         verbose = TRUE)

#saveRDS(fit,file="HRS_sparseVCBART.rds")

fit <- readRDS("HRS_sparseVCBART.rds")
assign(paste0("chain", chain_num),
       list(trees = fit$trees,
            varcounts = fit$varcounts[-(1:burn),,],
            sigma = fit$sigma))

save(list = paste0("chain", chain_num), file = paste0("hrs_chain", chain_num, ".RData"))

lam_draws <- fit$lambda 
lam_med   <- apply(lam_draws, 2, median)
lam_med_pred <- lam_med[-1]         
names(lam_med_pred) <- colnames(X1_all)

################## Re-fitting w/o INCOME & LBR ###########################
##########################################################################

chain_num <- 1
set.seed(424 + chain_num)

## ----- 1) Drop INCOME and LBR from predictors -----
drop_vars <- c("INCOME", "LBR")
stopifnot(all(drop_vars %in% colnames(X1_all)))

X1_reduced <- X1_all[, !(colnames(X1_all) %in% drop_vars), drop = FALSE]

p <- ncol(X1_reduced)

R_cont <- ncol(Z_cont_all)
R_cat  <- ncol(Z_cat_all)
R <- R_cont + R_cat

M  <- 50
mu0 <- rep(0, times = p + 1)
tau <- rep(0.5 / sqrt(M), times = p + 1)

## ----- 3) Build subject IDs (as before) -----
subj_id_all <- integer(0)
for (i in 1:n_all) {
  subj_id_all <- c(subj_id_all, rep(i, times = ni_all[i]))
}

## ----- 4) Fit sparseVCBART with reduced X -----
nd   <- 1000
burn <- 1000
thin <- 1

fit_reduced <- sparseVCBART::VCBART_cs(
  Y_train        = Y_all,
  subj_id_train  = subj_id_all,
  ni_train       = ni_all,
  X_train        = X1_reduced,          
  Z_cont_train   = Z_cont_all,
  Z_cat_train    = Z_cat_all,
  unif_cuts      = unif_cuts,
  cutpoints_list = cutpoints_list,
  cat_levels_list = cat_levels_list,
  sparse = TRUE,
  M = M, mu0 = mu0, tau = tau,
  nd = nd, burn = burn, thin = thin,
  save_trees = TRUE, save_samples = TRUE,
  verbose = TRUE
)

#saveRDS(fit_reduced, file = sprintf("HRS_sparseVCBART_noINCOME_noLBR.rds"))

assign(paste0("chain", chain_num),
       list(trees = fit_reduced$trees,
            varcounts = fit_reduced$varcounts[-(1:burn), , , drop = FALSE],
            sigma = fit_reduced$sigma))

save(list = paste0("chain", chain_num),
     file = sprintf("hrs_noINCOME_noLBR_chain%d.RData", chain_num))

lam_draws <- fit_reduced$lambda
lam_med   <- apply(lam_draws, 2, median)

lam_med_pred <- lam_med[-1]                 # drop intercept
names(lam_med_pred) <- colnames(X1_reduced)

############################# Lambda plots for full and reduced models ################################
#######################################################################################################
tidy_lambda <- function(fit, pred_names, scenario) {
  lam <- fit$lambda[, -1, drop = FALSE]   # drop intercept
  colnames(lam) <- pred_names
  
  as.data.frame(lam) %>%
    pivot_longer(everything(), names_to = "covariate", values_to = "lambda") %>%
    mutate(scenario = scenario)
}

df_full <- tidy_lambda(fit_full, colnames(X1_all),     "Full Model (p = 5)")
df_red  <- tidy_lambda(fit_red,  colnames(X1_reduced), "Reduced Model (p = 3)")
df <- bind_rows(df_full, df_red)

# consistent order (unused levels will drop within each facet if scales="free_x")
cov_order <- c("CHLD_SES","HS_COMPL","COLLEGE_COMPL","INCOME","LBR")
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
