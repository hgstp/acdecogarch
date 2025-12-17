
#' Generate irregular times with exponential gaps
#' @export
irr_times_exp <- function(lambda, t){ gaps <- stats::rexp(t, rate=lambda); c(0, cumsum(gaps)) }
#' Make matrix entries non-negative by absolute value
#' @export
pos_matrix <- function(A){ abs(A) }
#' Plot autocorrelation with 95% CI
#' @export
plot_acf_with_ci <- function(x, n){ ac <- stats::acf(x, plot=FALSE, lag.max=n); ci <- 1.96/sqrt(length(x));
  graphics::plot(0:n, c(1, ac$acf[-1]), type='h', xlab='Lag', ylab='ACF'); graphics::abline(h=0);
  graphics::abline(h=c(ci,-ci), lty=2, col='red') }
