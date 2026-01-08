#' Ljung-Box test on squared residuals for each column of a return matrix
#' @export
lb_residuals_matrix <- function(G, sigma, alpha = 0.05, lag = NULL) {
  stopifnot(is.matrix(G), is.matrix(sigma))
  stopifnot(nrow(G) == nrow(sigma))
  n <- nrow(G)
  p <- ncol(G)
  if (is.null(lag)) lag <- max(1L, floor(sqrt(n)))
  rows <- lapply(seq_len(p), function(k) {
    res <- G[, k] / sqrt(sigma[, k])
    ht <- stats::Box.test(res^2, lag = lag, type = "Ljung-Box")
    data.frame(series = k, H = as.integer(ht$p.value < alpha), p.value = ht$p.value, Q = ht$statistic)
  })
  do.call(rbind, rows)
}
