summarize_beta <- function(beta_samples,
                           level = 0.95,
                           include_median = FALSE,
                           na.rm = FALSE) {
  dims <- dim(beta_samples)
  if (length(dims) != 3L) {
    stop("beta_samples must be a 3D array of shape nd x n x p.")
  }
  nd <- dims[1]; n <- dims[2]; p <- dims[3]

  if (level <= 0 || level >= 1) stop("level must be in (0,1).")
  alpha  <- (1 - level) / 2
  probs  <- c(alpha, 1 - alpha)
  qnames <- paste0(c("L", "U"), round(level * 100))

  stat_names <- if (include_median) c("MEAN", qnames[1], "MEDIAN", qnames[2]) else c("MEAN", qnames[1], qnames[2])
  out <- array(NA_real_, dim = c(n, length(stat_names), p),
               dimnames = list(NULL, stat_names, NULL))

  for (j in seq_len(p)) {
    tmp <- beta_samples[ , , j, drop = FALSE]
    tmp <- matrix(tmp, nrow = nd, ncol = n)  

    out[, "MEAN", j] <- colMeans(tmp, na.rm = na.rm)

    q <- apply(tmp, 2L, stats::quantile, probs = probs, na.rm = na.rm, names = FALSE)
    out[, qnames[1], j] <- q[1, ]
    out[, qnames[2], j] <- q[2, ]

    if (include_median) {
      out[, "MEDIAN", j] <- apply(tmp, 2L, stats::median, na.rm = na.rm)
    }
  }
  out
}
