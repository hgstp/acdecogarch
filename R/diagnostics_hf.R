#' Ljung-Box test for HF residuals
#' @export
lb_hf <- function(ret, sigma, alpha = 0.05, lag = NULL) {
  stopifnot(length(ret) == length(sigma))
  res <- ret / sqrt(sigma)
  n <- length(res)
  if (is.null(lag)) lag <- max(1L, floor(sqrt(n)))
  ht <- stats::Box.test(res^2, lag = lag, type = "Ljung-Box")
  crit <- stats::qchisq(1 - alpha, df = lag)
  list(H = as.integer(ht$statistic > crit), p.value = ht$p.value, Q = ht$statistic, critical = crit)
}
