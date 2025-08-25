predict_betas <- function(fit,
                          Z_cont = matrix(0, nrow = 1, ncol = 1),
                          Z_cat  = matrix(0L, nrow = 1, ncol = 1),
                          verbose = TRUE) {

  if (is.null(fit$trees)) {
    stop("`fit$trees` is NULL. Refit with `save_trees = TRUE` to enable prediction.")
  }
  if (is.null(fit$x_mean) || is.null(fit$x_sd) ||
      is.null(fit$y_mean) || is.null(fit$y_sd)) {
    stop("`fit` must contain x_mean, x_sd, y_mean, y_sd (produced by VCBART_*).")
  }

  if (!is.matrix(Z_cont)) stop("Z_cont must be a matrix.")
  if (!is.matrix(Z_cat))  stop("Z_cat must be a matrix.")
  if (!is.integer(Z_cat)) storage.mode(Z_cat) <- "integer"  # coerce if needed

  p <- length(fit$x_mean)

  M <- length(fit$trees[[1]][[1]])

  beta_draws_std <- .predict_vcbart(
    tree_draws = fit$trees,
    p = p, M = M,
    tZ_cont = t(Z_cont),
    tZ_cat  = t(Z_cat),
    verbose = verbose
  )

  beta_draws <- rescale_beta(
    beta_draws_std,
    fit$y_mean, fit$y_sd,
    fit$x_mean, fit$x_sd
  )

  beta_draws
}
