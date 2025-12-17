
#' Simulate univariate ECOGARCH(1,1) driven by compound Poisson
#' @export
simulate_ecogarch11_cp_univariate <- function(par, lambda=1, nmax=3000, nu=6, equidistant_step=NULL){
  a1 <- par$a1; theta <- par$theta; gamma <- par$gamma; mu <- par$mu
  dt <- stats::rexp(nmax, rate=lambda)
  t_jump <- cumsum(dt)
  Z <- stats::rt(nmax, df=nu) * sqrt((nu-2)/(nu*lambda))
  K <- (sqrt((nu-2)/pi) * gamma(0.5*(nu-1))/gamma(0.5*nu)) / sqrt(lambda)
  X <- numeric(nmax); s2_left <- numeric(nmax)
  for (i in seq_len(nmax)){
    Xprev <- if (i==1) 0 else X[i-1]; dti <- dt[i]
    s2_left[i] <- exp(mu + exp(-a1*dti)*Xprev - gamma*lambda*K*(1 - exp(-a1*dti))/abs(a1))
    X[i] <- exp(-a1*dti)*Xprev + theta*Z[i] + gamma*(abs(Z[i]) - lambda*K*(1 - exp(-a1*dti))/abs(a1))
  }
  G <- cumsum(sqrt(s2_left) * Z)
  if (is.null(equidistant_step)) return(list(G=G, t_obs=t_jump, vol=exp(mu + X)))
  t_grid <- seq(equidistant_step, max(t_jump), by=equidistant_step)
  all_times <- sorted <- sort(c(t_jump, t_grid))
  Z_all <- rep(0, length(all_times)); Z_all[match(t_jump, all_times)] <- Z
  s2_all <- numeric(length(all_times)); X_all <- numeric(length(all_times))
  last_t <- 0; Xcur <- 0
  for (i in seq_along(all_times)){
    dti <- all_times[i] - last_t
    s2_all[i] <- exp(mu + exp(-a1*dti)*Xcur - gamma*lambda*K*(1 - exp(-a1*dti))/abs(a1))
    Zi <- Z_all[i]
    Xcur <- exp(-a1*dti)*Xcur + theta*Zi + gamma*(abs(Zi) - lambda*K*(1 - exp(-a1*dti))/abs(a1))
    X_all[i] <- Xcur; last_t <- all_times[i]
  }
  G_all <- cumsum(sqrt(s2_all) * Z_all)
  list(G=G_all, t_obs=all_times, vol=exp(mu + X_all))
}
