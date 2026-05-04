############################ Yeast cell-cycle data analysis #############################
#########################################################################################
library(sparseVCBART)
library(dplyr)
yeast <- read.csv("~/sparseVCBART/yeast.csv")
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
make_cat_levels <- function(Z_cat) {
  if (!is.matrix(Z_cat)) Z_cat <- as.matrix(Z_cat)
  lapply(seq_len(ncol(Z_cat)), function(j) sort(unique(Z_cat[, j])))
}

stopifnot(all(c("mRNA","Time","Gene") %in% names(yeast)))
tf_cols <- setdiff(names(yeast), c("mRNA","Time","Gene"))

yeast <- yeast[order(yeast$Gene, yeast$Time), ]

N <- nrow(yeast)
cat("N =", N,
    " n_genes =", length(unique(yeast$Gene)),
    " n_times =", length(unique(yeast$Time)), "\n")

subj_id_all <- as.integer(factor(yeast$Gene))  
r <- rle(subj_id_all)
stopifnot(all(r$values == seq_len(max(subj_id_all))))
ni_all <- as.integer(r$lengths)

#Y_all <- as.numeric(scale(yeast$mRNA)) 
Y_all <- as.numeric(yeast$mRNA) 

X1_all <- as.matrix(yeast[, tf_cols, drop = FALSE])

#X1_all <- scale(X1_all)

t_raw <- yeast$Time
t_scaled <- 2 * (t_raw - min(t_raw)) / (max(t_raw) - min(t_raw)) - 1

## Z_cont: must be a numeric matrix
Z_cont_all <- matrix(as.numeric(t_scaled), ncol = 1)
colnames(Z_cont_all) <- "Time"

## No categorical modifiers: use NULL (NOT a 0-col matrix)
Z_cat_all <- matrix(0, nrow = 1, ncol = 1)
cp_list <- make_cutpoints(Z_cont_all, n_cuts = 50)
cat_levels_list <- NULL   # or list()

stopifnot(length(Y_all) == nrow(X1_all),
          length(Y_all) == nrow(Z_cont_all))

fit_sparse_yeast <- sparseVCBART::VCBART_cs(
  Y_train         = Y_all,
  subj_id_train   = subj_id_all,
  ni_train        = ni_all,
  X_train         = X1_all,
  Z_cont_train    = Z_cont_all,
  Z_cat_train     = Z_cat_all,       # Passed as 1x1 zero matrix
  cutpoints_list  = cp_list,
  cat_levels_list = cat_levels_list, # Passed as NULL
  M = 50, nd = 2000, burn = 500,
  save_samples = TRUE, save_trees = TRUE, verbose = TRUE
)

lam_draws <- fit_sparse_yeast$lambda   # [iter x (p+1)] includes intercept
lam_med   <- apply(lam_draws, 2, median)

lam_med_pred <- lam_med[-1]
names(lam_med_pred) <- colnames(X1_all)

sort(lam_med_pred, decreasing = TRUE)[1:20]

saveRDS(fit_sparse_yeast,file="fit_sparseVCBART_yeast.rds")
fit_sparse_yeast <- readRDS("fit_sparseVCBART_yeast.rds")

####################### Junk modifiers ######################################
t_raw <- yeast$Time
t_scaled <- 2 * (t_raw - min(t_raw)) / (max(t_raw) - min(t_raw)) - 1

# Generate 20 spurious modifiers drawn from U(-1, 1)
set.seed(123) # Essential for reproducibility of the robustness check
n_spurious <- 19
Z_spurious <- matrix(runif(N * n_spurious, min = -1, max = 1), 
                     nrow = N, ncol = n_spurious)
colnames(Z_spurious) <- paste0("Junk_", 1:n_spurious)

## Z_cont: must be a numeric matrix (Now 21 columns)
Z_cont_all <- cbind(Time = as.numeric(t_scaled), Z_spurious)

## No categorical modifiers
Z_cat_all <- matrix(0, nrow = 1, ncol = 1)
cp_list <- make_cutpoints(Z_cont_all, n_cuts = 50)
cat_levels_list <- NULL   

stopifnot(length(Y_all) == nrow(X1_all),
          length(Y_all) == nrow(Z_cont_all))

# Fit the robust model
fit_sparse_yeast_robust <- sparseVCBART::VCBART_cs(
  Y_train         = Y_all,
  subj_id_train   = subj_id_all,
  ni_train        = ni_all, # Note: using your ni_all variable
  X_train         = X1_all,
  Z_cont_train    = Z_cont_all,
  Z_cat_train     = Z_cat_all,       
  cutpoints_list  = cp_list,
  cat_levels_list = cat_levels_list, 
  M = 30, nd = 2000, burn = 500,
  sparse = TRUE,
  save_samples = TRUE, save_trees = TRUE, verbose = TRUE
)

saveRDS(fit_sparse_yeast_robust,file="fit_sparseVCBART_yeast_augmented.rds")

theta_full <- fit_sparse_yeast_robust$theta

# 2. Discard the first 500 burn-in iterations 
# Dimensions change from [2500, 20, 97] to [2000, 20, 97]
theta_post_burn <- theta_full[501:2500, , ]

# 3. Calculate the posterior Mean across iterations (Dimension 1)
# MARGIN = c(2, 3) tells R to keep the Modifier and Covariate dimensions intact
theta_mean_matrix <- apply(theta_post_burn, MARGIN = c(2, 3), FUN = mean)

# 4. Calculate the posterior Median across iterations (Dimension 1)
theta_median_matrix <- apply(theta_post_burn, MARGIN = c(2, 3), FUN = median)

# idx_time <- 2
# theta_time_vector <- theta_median_matrix[idx_time, ]
# 
# # 3. Summarize this vector to a single global value
# # The mean gives the average modifier importance across all transcription factors
# time_global_mean <- mean(theta_time_vector)
# 
# # Alternatively, the median gives the strictly middle importance
# time_global_median <- median(theta_time_vector)
# 
# # You can also look at the maximum to see its strongest single effect
# time_global_max <- max(theta_time_vector)
# 
# cat("Global Mean Theta for Time:  ", time_global_mean, "\n")
# cat("Global Median Theta for Time:", time_global_median, "\n")
# cat("Maximum Theta for Time:      ", time_global_max, "\n")
# 
# 
# # 1. Identify the rows for the spurious modifiers (everything EXCEPT Time)
# # Assuming idx_time is the row index for Time as defined previously
# idx_junk <- setdiff(seq_len(nrow(theta_median_matrix)), idx_time)[-1]
# 
# # 2. Extract the sub-matrix containing ONLY the junk modifiers
# theta_junk_matrix <- theta_median_matrix[idx_junk, ]
# junk_global_means <- rowMeans(theta_junk_matrix)
# max_junk_global_mean <- max(junk_global_means)
# top_junk_name <- rownames(theta_junk_matrix)[which.max(junk_global_means)]
# 
# cat("Global Mean Theta for Time:             ", time_global_mean, "\n")
# cat("Maximum Global Mean among Junk Modifiers:", max_junk_global_mean, "\n")
# cat("Top Junk Modifier:                      ", top_junk_name, "\n")
# 
# pct_diff <- ((time_global_mean - max_junk_global_mean) / max_junk_global_mean) * 100
# cat("Time is higher by:                      ", round(pct_diff, 2), "%\n")


# 1. Assign names so we never have to guess indices again
rownames(theta_median_matrix) <- colnames(Z_cont_all)

# 2. Safely extract Time by name
idx_time <- which(rownames(theta_median_matrix) == "Time")
theta_time_vector <- theta_median_matrix[idx_time, ]
time_global_mean <- mean(theta_time_vector)

# 3. Safely extract all Junk variables by name
idx_junk <- which(rownames(theta_median_matrix) != "Time")
theta_junk_matrix <- theta_median_matrix[idx_junk, ]
junk_global_means <- rowMeans(theta_junk_matrix)

# 4. Find the absolute maximum global mean among the junk variables
max_junk_global_mean <- max(junk_global_means)
top_junk_name <- rownames(theta_junk_matrix)[which.max(junk_global_means)]

# Print the corrected results
cat("Global Mean Theta for Time:              ", time_global_mean, "\n")
cat("Maximum Global Mean among Junk Modifiers:", max_junk_global_mean, "\n")
cat("Top Junk Modifier:                       ", top_junk_name, "\n")

# Calculate the actual percentage difference
pct_diff <- ((time_global_mean - max_junk_global_mean) / max_junk_global_mean) * 100
cat("Time is higher by:                       ", round(pct_diff, 2), "%\n")

######################### varcounts #######################################
vc_post_burn <- fit_sparse_yeast_robust$varcounts[501:2500, , ]
vc_summary_matrix <- apply(vc_post_burn, MARGIN = c(2, 3), FUN = mean)

# Now this will work perfectly (20 names for 20 rows)
rownames(vc_summary_matrix) <- colnames(Z_cont_all)

# 2. Safely extract Time by name
idx_time <- which(rownames(vc_summary_matrix) == "Time")
vc_time_vector <- vc_summary_matrix[idx_time, ]

# Summarize to a single global value (Mean splits per covariate)
time_vc_global <- mean(vc_time_vector) 

# 3. Safely extract all Junk variables by name
idx_junk <- which(rownames(vc_summary_matrix) != "Time")
vc_junk_matrix <- vc_summary_matrix[idx_junk, ]

# Calculate the global mean varcount for each junk modifier
junk_vc_global <- rowMeans(vc_junk_matrix) 

# 4. Find the absolute maximum global varcount among the junk variables
max_junk_vc_global <- max(junk_vc_global)
top_junk_vc_name <- rownames(vc_junk_matrix)[which.max(junk_vc_global)]

# Print the results
cat("Global Mean VarCount for Time:             ", time_vc_global, "\n")
cat("Maximum Global Mean among Junk Modifiers:  ", max_junk_vc_global, "\n")
cat("Top Junk Modifier (VarCounts):             ", top_junk_vc_name, "\n")

# Calculate the percentage difference
pct_diff_vc <- ((time_vc_global - max_junk_vc_global) / max_junk_vc_global) * 100
cat("Time VarCounts are higher by:              ", round(pct_diff_vc, 2), "%\n")
################################################################################
time_grid_raw <- sort(unique(yeast$Time))
tmin <- min(time_grid_raw); tmax <- max(time_grid_raw)
time_grid_scaled <- 2*(time_grid_raw - tmin)/(tmax - tmin) - 1

Z_cont_grid <- matrix(time_grid_scaled, ncol = 1)
colnames(Z_cont_grid) <- "Time"

Z_cat_grid <- Z_cat_all[rep(1, nrow(Z_cont_grid)), , drop = FALSE]

# b_draws <- sparseVCBART::predict_betas(
#   fit_sparse_yeast,
#   Z_cont = Z_cont_grid,
#   Z_cat  = Z_cat_grid,
#   verbose = TRUE
# )
# b_sum <- sparseVCBART::summarize_beta(b_draws)  # [grid x {MEAN,L95,U95} x (p+1)]
# topJ <- order(lam_med_pred, decreasing=TRUE)[1:3]  # top 6 TFs
# par(mfrow=c(3,2))
# for (j in topJ) {
#   mu  <- b_sum[, "MEAN", j+1]  # +1 because intercept is column 1
#   lo  <- b_sum[, "L95",  j+1]
#   hi  <- b_sum[, "U95",  j+1]
#   plot(time_grid_raw, mu, type="l",
#        main=colnames(X1_all)[j], xlab="Time", ylab=expression(beta(t)))
#   lines(time_grid_raw, lo, lty=2); lines(time_grid_raw, hi, lty=2)
#   abline(h=0, lty=3)
# }
# par(mfrow=c(1,1))


# 1. Get predictions and summarize
b_draws <- sparseVCBART::predict_betas(
  fit_sparse_yeast,
  Z_cont = Z_cont_grid,
  Z_cat  = Z_cat_grid,
  verbose = TRUE
)
b_sum <- sparseVCBART::summarize_beta(b_draws)  # [grid x {MEAN,L95,U95} x (p+1)]

# 2. Identify the top 3 TFs
#topJ <- order(lam_med_pred, decreasing = TRUE)[1:3]

topJ <- order(lam_med_pred, decreasing = TRUE)[4:7]

# 3. Restructure the data into a single data frame for ggplot
plot_data <- do.call(rbind, lapply(topJ, function(j) {
  data.frame(
    Time  = time_grid_raw,
    Mean  = b_sum[, "MEAN", j + 1],
    Lower = b_sum[, "L95",  j + 1],
    Upper = b_sum[, "U95",  j + 1],
    TF    = colnames(X1_all)[j]
  )
}))

# Convert TF to a factor so ggplot keeps them ordered by importance (1st, 2nd, 3rd)
plot_data$TF <- factor(plot_data$TF, levels = colnames(X1_all)[topJ])

# 4. Plot using ggplot2
p <- ggplot(plot_data, aes(x = Time)) +
  # Add horizontal baseline at 0
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  
  # Optional: Shaded region for credible intervals 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.2) +
  
  # Bolder mean trajectory lines
  geom_line(aes(y = Mean), color = "steelblue", linewidth = 1.2) +
  
  # Dashed lines for the credible intervals (matching your base R approach)
  geom_line(aes(y = Lower), color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = Upper), color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  
  # Create a 1x3 grid of plots
  facet_wrap(~ TF, ncol = 4, scales = "free_y") +
  
  # Labels and themes
  labs(x = "Time", y = expression(beta(t))) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(p)

ggsave("TEs_yeast_other.png", plot = p, width = 11, height = 4, dpi = 300)


############################### Running 4 chains ############################
# --- settings ---
M    <- 50
nd   <- 2000
burn <- 400
thin <- 1

n_chains  <- 4
base_seed <- 1001
seeds <- base_seed + 1:n_chains

fit_one_chain <- function(seed) {
  set.seed(seed)
  sparseVCBART::VCBART_cs(
    Y_train         = Y_all,
    subj_id_train   = subj_id_all,
    ni_train        = ni_all,
    X_train         = X1_all,
    Z_cont_train    = Z_cont_all,
    Z_cat_train     = Z_cat_all,        # your 1x1 dummy
    cutpoints_list  = cp_list,
    cat_levels_list = cat_levels_list,  # NULL
    M = M, nd = nd, burn = burn, thin = thin,
    save_samples = TRUE,
    save_trees   = TRUE,               # recommended for multi-chain
    verbose      = TRUE
  )
}

fits <- vector("list", n_chains)
for (cc in 1:n_chains) {
  cat("\n--- Fitting chain", cc, "seed", seeds[cc], "---\n")
  fits[[cc]] <- fit_one_chain(seeds[cc])
}

# (optional) save fits
saveRDS(fits, "yeast_sparseVCBART_4chains.rds")

# lambda: [iter x (p+1)] (includes intercept in column 1)
lam_list <- lapply(fits, function(f) f$lambda)
lam_pool <- do.call(rbind, lam_list)  # pooled draws across chains

lam_med <- apply(lam_pool, 2, median, na.rm = TRUE)
lam_med_pred <- lam_med[-1]  # drop intercept
names(lam_med_pred) <- colnames(X1_all)


# posterior medians
lam_med <- apply(lam_pool, 2, median, na.rm = TRUE)
lam_med_pred <- lam_med[-1]   # drop intercept
names(lam_med_pred) <- colnames(X1_all)

# build data frame and keep top 20
plot_df <- data.frame(
  Feature = names(lam_med_pred),
  MedianLambda = as.numeric(lam_med_pred)
) %>%
  arrange(desc(MedianLambda)) %>%
  slice_head(n = 20) %>%
  mutate(
    Group = ifelse(MedianLambda > 10, "Signal", "Noise"),
    Feature = factor(Feature, levels = Feature)
  )


# lollipop plot
yeast_imp <- ggplot(plot_df, aes(x = Feature, y = MedianLambda, color = Group)) +
  geom_segment(aes(xend = Feature, y = 0, yend = MedianLambda),
               linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Signal" = "forestgreen", "Noise" = "red3")) +
  labs(
    x = NULL,
    y = expression("median " * lambda[j]),
    color = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor.x = element_blank()
  )
ggsave("TEs_imp.png", plot = yeast_imp, width = 11, height = 5, dpi = 300)

# top TFs by median lambda
sort(lam_med_pred, decreasing = TRUE)[1:20]

lam_mean <- apply(lam_pool, 2, mean, na.rm = TRUE)
lam_mean_pred <- lam_mean[-1]  # drop intercept
names(lam_mean_pred) <- colnames(X1_all)

# top TFs by median lambda
sort(lam_mean_pred, decreasing = TRUE)[1:20]


time_grid_raw <- sort(unique(yeast$Time))
tmin <- min(time_grid_raw); tmax <- max(time_grid_raw)
time_grid_scaled <- 2*(time_grid_raw - tmin)/(tmax - tmin) - 1

Z_cont_grid <- matrix(time_grid_scaled, ncol = 1)
colnames(Z_cont_grid) <- "Time"

# replicate the 1x1 dummy categorical matrix to grid length
Z_cat_grid <- Z_cat_all[rep(1, nrow(Z_cont_grid)), , drop = FALSE]

if (!requireNamespace("abind", quietly = TRUE)) install.packages("abind")

b_draws_list <- lapply(fits, function(f) {
  sparseVCBART::predict_betas(
    f,
    Z_cont = Z_cont_grid,
    Z_cat  = Z_cat_grid,
    verbose = TRUE
  )
})

# Pool draws by stacking iterations across chains
# (assumes first dim is iteration)
b_draws_pool <- abind::abind(b_draws_list, along = 1)

# Summarize pooled draws: [grid x {MEAN,L95,U95} x (p+1)]
b_sum_pool <- sparseVCBART::summarize_beta(b_draws_pool)

# topJ <- order(lam_med_pred, decreasing = TRUE)[1:3]  # top 6 TFs
# 
# par(mfrow = c(1,3))
# for (j in topJ) {
#   mu <- b_sum_pool[, "MEAN", j + 1]  # +1 for intercept
#   lo <- b_sum_pool[, "L95",  j + 1]
#   hi <- b_sum_pool[, "U95",  j + 1]
#   
#   plot(time_grid_raw, mu, type = "l",
#        main = colnames(X1_all)[j],
#        xlab = "Time", ylab = expression(beta(t)))
#   lines(time_grid_raw, lo, lty = 2)
#   lines(time_grid_raw, hi, lty = 2)
#   abline(h = 0, lty = 3)
# }
# par(mfrow = c(1,1))


topJ <- order(lam_med_pred, decreasing = TRUE)[4:10]

# 3. Restructure the data into a single data frame for ggplot
plot_data <- do.call(rbind, lapply(topJ, function(j) {
  data.frame(
    Time  = time_grid_raw,
    Mean  = b_sum_pool[, "MEAN", j + 1],
    Lower = b_sum_pool[, "L95",  j + 1],
    Upper = b_sum_pool[, "U95",  j + 1],
    TF    = colnames(X1_all)[j]
  )
}))

# Convert TF to a factor so ggplot keeps them ordered by importance (1st, 2nd, 3rd)
plot_data$TF <- factor(plot_data$TF, levels = colnames(X1_all)[topJ])

# 4. Plot using ggplot2
p <- ggplot(plot_data, aes(x = Time)) +
  # Add horizontal baseline at 0
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8) +
  
  # Optional: Shaded region for credible intervals 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.2) +
  
  # Bolder mean trajectory lines
  geom_line(aes(y = Mean), color = "steelblue", linewidth = 1.2) +
  
  # Dashed lines for the credible intervals (matching your base R approach)
  geom_line(aes(y = Lower), color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = Upper), color = "steelblue", linetype = "dashed", linewidth = 0.8) +
  
  # Create a 1x3 grid of plots
  facet_wrap(~ TF,nrow=2, ncol = 4, scales = "free_y") +
  
  # Labels and themes
  labs(x = "Time", y = expression(beta(t))) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

print(p)

#ggsave("TEs_yeast_imp.png", plot = p, width = 11, height = 4, dpi = 300)
ggsave("TEs_yeast_rem.png", plot = p, width = 11, height = 5, dpi = 300)

#################################################################################################
################################ Comparing MSPEs ################################################

## Combine posterior mean betas into posterior mean of E[Y | X, Z]
yhat_from_beta_summ <- function(sum_obj, X_test) {
  beta0 <- sum_obj[, "MEAN", 1]
  if (ncol(X_test) == 0) return(as.numeric(beta0))
  B <- sapply(1:ncol(X_test), function(j) sum_obj[, "MEAN", j + 1])  # [n_test x p]
  as.numeric(beta0 + rowSums(X_test * B))
}

build_vcbart_inputs <- function(df, tf_cols, tmin, tmax,
                                y_mu = NULL, y_sd = NULL,
                                x_mu = NULL, x_sd = NULL) {
  df <- df[order(df$Gene, df$Time), ]
  
  ## subject ids must be consecutive within each split
  subj_id <- as.integer(factor(df$Gene))
  ni      <- as.integer(rle(subj_id)$lengths)
  
  ## Y scaling (use train stats if provided)
  y_raw <- df$mRNA
  if (is.null(y_mu)) y_mu <- mean(y_raw)
  if (is.null(y_sd)) y_sd <- sd(y_raw)
  #y <- (y_raw - y_mu) / y_sd
  y <- y_raw
  ## X scaling (use train stats if provided)
  X_raw <- as.matrix(df[, tf_cols, drop = FALSE])
  # if (is.null(x_mu)) x_mu <- colMeans(X_raw)
  # if (is.null(x_sd)) x_sd <- apply(X_raw, 2, sd)
  # stopifnot(all(x_sd > 0))
  # X <- sweep(sweep(X_raw, 2, x_mu, "-"), 2, x_sd, "/")
  X <- X_raw
  ## Z_cont: scaled time in [-1,1] (use global min/max since grid is fixed anyway)
  t_scaled <- 2 * (df$Time - tmin) / (tmax - tmin) - 1
  Z_cont <- matrix(as.numeric(t_scaled), ncol = 1)
  colnames(Z_cont) <- "Time"
  
  list(
    Y = as.numeric(y),
    X = X,
    Z_cont = Z_cont,
    subj_id = subj_id,
    ni = ni,
    y_mu = y_mu, y_sd = y_sd,
    x_mu = x_mu, x_sd = x_sd
  )
}
one_split_mspe <- function(yeast, tf_cols, Z_cat_dummy,
                           M = 50, nd = 1200, burn = 400,
                           n_train_genes = 37) {
  
  genes <- sort(unique(yeast$Gene))
  train_genes <- sample(genes, n_train_genes)
  test_genes  <- setdiff(genes, train_genes)
  
  df_tr <- yeast[yeast$Gene %in% train_genes, ]
  df_te <- yeast[yeast$Gene %in% test_genes,  ]
  
  tgrid <- sort(unique(yeast$Time))
  tmin <- min(tgrid); tmax <- max(tgrid)
  
  tr <- build_vcbart_inputs(df_tr, tf_cols, tmin, tmax)
  te <- build_vcbart_inputs(df_te, tf_cols, tmin, tmax,
                            y_mu = tr$y_mu, y_sd = tr$y_sd,
                            x_mu = tr$x_mu, x_sd = tr$x_sd)
  
  ## cutpoints from training Z
  cp_list <- make_cutpoints(tr$Z_cont, n_cuts = 50)
  
  fit <- sparseVCBART::VCBART_cs(
    Y_train         = tr$Y,
    subj_id_train   = tr$subj_id,
    ni_train        = tr$ni,
    X_train         = tr$X,
    Z_cont_train    = tr$Z_cont,
    Z_cat_train     = Z_cat_dummy,     # your 1x1 zero matrix
    cutpoints_list  = cp_list,
    cat_levels_list = NULL,
    M = M, nd = nd, burn = burn,
    save_samples = TRUE,
    save_trees   = TRUE,
    verbose = TRUE
  )
  
  ## predict betas at each test-row's (time) modifier
  Z_cat_te <- Z_cat_dummy[rep(1, nrow(te$Z_cont)), , drop = FALSE]
  
  b_te <- sparseVCBART::predict_betas(
    fit, Z_cont = te$Z_cont, Z_cat = Z_cat_te, verbose = FALSE
  )
  sum_te <- sparseVCBART::summarize_beta(b_te)
  
  yhat <- yhat_from_beta_summ(sum_te, te$X)
  
  ## MSPE on standardized scale (matches how we fit/predict)
  mspe <- mean((yhat - te$Y)^2)
  
  ## baseline: intercept-only (time-varying intercept, no TFs)
  mspe0 <- mean((sum_te[, "MEAN", 1] - te$Y)^2)
  
  c(mspe = mspe, mspe_intercept_only = mspe0)
}

set.seed(1)

R <- 20
out <- matrix(NA_real_, nrow = R, ncol = 2,
              dimnames = list(NULL, c("mspe","mspe_intercept_only")))

for (r in 1:R) {
  out[r, ] <- one_split_mspe(
    yeast = yeast,
    tf_cols = tf_cols,
    Z_cat_dummy = Z_cat_all,   # your 1x1 zero matrix
    M = 50, nd = 1200, burn = 400,
    n_train_genes = 37
  )
  if (r %% 10 == 0) cat("done", r, "of", R, "\n")
}

mspe_mean <- mean(out[, "mspe"])
mspe_se   <- sd(out[, "mspe"]) / sqrt(R)

base_mean <- mean(out[, "mspe_intercept_only"])
base_se   <- sd(out[, "mspe_intercept_only"]) / sqrt(R)

cat("sparseVCBART MSPE:", mspe_mean, " (SE:", mspe_se, ")\n")
cat("Intercept-only MSPE:", base_mean, " (SE:", base_se, ")\n")

