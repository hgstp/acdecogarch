#' Estimate volatility and innovations given ECOGARCH(1,1) parameters
#' @param par list(a1, theta, gamma, mu, lambda, K)
#' @param G log-prices
#' @param tQML observation times
#' @param method 'approx'|'exact'
#' @export
estimate_sigma_ecogarch <- function(par, G, tQML, method = c("approx", "exact")) {
  method <- match.arg(method)
  R <- c(G[1], diff(G))
  dt <- c(tQML[1], diff(tQML))
  n <- length(R)
  X <- numeric(n)
  s2 <- numeric(n)
  Zhat <- numeric(n)
  for (i in seq_len(n)) {
    Xprev <- if (i == 1) 0 else X[i - 1]
    if (method == "approx") {
      s2[i] <- exp(par$mu + exp(-par$a1 * dt[i]) * Xprev - par$gamma * par$lambda * par$K * dt[i])
      Zhat[i] <- R[i] / sqrt(s2[i])
      X[i] <- exp(-par$a1 * dt[i]) * Xprev + par$theta * Zhat[i] + par$gamma * abs(Zhat[i]) - par$gamma * par$lambda * par$K * dt[i]
    } else {
      s2[i] <- exp(par$mu + exp(-par$a1 * dt[i]) * Xprev - par$gamma * par$lambda * par$K * (1 - exp(-par$a1 * dt[i])) / par$a1)
      Zhat[i] <- R[i] / sqrt(s2[i])
      X[i] <- exp(-par$a1 * dt[i]) * Xprev + par$theta * Zhat[i] + par$gamma * (abs(Zhat[i]) - par$lambda * par$K * (1 - exp(-par$a1 * dt[i])) / par$a1)
    }
  }
  list(sigma = s2, innovations = Zhat)
}
