
#' Stepwise one-step prediction and coverage
#' @export
prediction_stepwise <- function(tQML, G, alpha=0.05, start_k=100, method='approx', grad='none'){
  R <- c(G[1], diff(G)); dt <- c(tQML[1], diff(tQML)); N <- length(R)
  sigma_hat <- rep(NA_real_, N); G_upper <- rep(NA_real_, N); G_lower <- rep(NA_real_, N)
  z <- stats::qnorm(1 - alpha/2)
  for (k in seq.int(start_k, N-1)){
    fit <- ecogarch_qmle(Gdt=R[1:k], dt=dt[1:k], init=list(a1=0.1, theta=-0.1, gamma=0.3, mu=-3), method=method, grad=grad)
    est <- estimate_sigma_ecogarch(fit$par, G[1:k+1], tQML[1:k+1], method)
    s2_next <- est$sigma[k+1]; sigma_hat[k+1] <- s2_next
    G_upper[k+1] <- G[k+1] + z * sqrt(s2_next / fit$par$lambda)
    G_lower[k+1] <- G[k+1] - z * sqrt(s2_next / fit$par$lambda)
  }
  idx <- which(!is.na(G_upper)); L <- length(idx)
  thirds <- max(1L, floor(L/3))
  seg1 <- idx[1:thirds]
  seg2 <- idx[(thirds+1):min(2*thirds, L)]
  seg3 <- idx[(min(2*thirds, L)+1):L]
  inside <- function(seg){ if(length(seg)) sum(G[seg] <= G_upper[seg] & G[seg] >= G_lower[seg]) else NA_integer_ }
  coverage <- c(inside(seg1), inside(seg2), inside(seg3))
  pct <- c(if(length(seg1)) coverage[1]/length(seg1) else NA_real_,
           if(length(seg2)) coverage[2]/length(seg2) else NA_real_,
           if(length(seg3)) coverage[3]/length(seg3) else NA_real_)
  list(G_upper=G_upper, G_lower=G_lower, sigma_hat=sigma_hat, coverage=coverage, coverage_pct=pct)
}

#' Coverage helper
#' @export
prediction_interval_coverage <- function(G, G_upper, G_lower, start_idx){
  idx <- seq.int(start_idx, length(G)); L <- length(idx); thirds <- floor(L/3)
  seg1 <- idx[1:thirds]; seg2 <- idx[(thirds+1):(2*thirds)]; seg3 <- idx[(2*thirds+1):L]
  inside <- function(seg){ if(length(seg)) sum(G[seg] <= G_upper[seg] & G[seg] >= G_lower[seg]) else NA_integer_ }
  counts <- c(inside(seg1), inside(seg2), inside(seg3))
  percents <- counts / c(length(seg1), length(seg2), length(seg3))
  list(counts=counts, percents=percents)
}
