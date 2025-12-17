
#' Previous-tick interpolation
#' @export
prev_tick <- function(t, x, X){ idx <- findInterval(t, x, left.open=FALSE); X[pmax(idx,1)] }

#' Piecewise lognormal one-step volatility density
#' @export
vol_prediction_density_piecewise <- function(Xn, dt_next, par, grid=NULL, exact=TRUE){
  a1 <- par$a1; theta <- par$theta; gamma <- par$gamma; mu <- par$mu; lam <- par$lambda; K <- par$K
  if (is.null(grid)) grid <- exp(seq(log(1e-6), log(100), length.out=1024))
  comp <- if (exact) (1 - exp(-a1*dt_next))/abs(a1) else dt_next
  varphi <- exp(-a1*dt_next)*Xn - gamma*lam*K*comp
  z_m <- (log(grid) - mu - varphi)/(theta + gamma); z_n <- (log(grid) - mu - varphi)/(theta - gamma)
  fZ <- function(z) sqrt(lam/(2*pi)) * exp(-0.5*lam*z^2)
  dens_m <- (1/grid) * fZ(z_m) * (1/abs(theta + gamma)); dens_n <- (1/grid) * fZ(z_n) * (1/abs(theta - gamma))
  m_idx <- which(log(grid) - mu < varphi); n_idx <- which(log(grid) - mu >= varphi)
  dens <- numeric(length(grid))
  if (theta < -gamma && -gamma < 0){ dens[m_idx] <- dens_m[m_idx]; dens[n_idx] <- dens_n[n_idx]
  } else if (-gamma < theta && theta < 0){ dens[m_idx] <- 0; dens[n_idx] <- dens_m[n_idx] + dens_n[n_idx]
  } else if (0 < gamma && gamma < theta){ dens[m_idx] <- dens_n[m_idx]; dens[n_idx] <- dens_m[n_idx]
  } else if (0 < theta && theta < gamma){ dens[m_idx] <- 0; dens[n_idx] <- dens_m[n_idx] + dens_n[n_idx]
  } else { dens[m_idx] <- dens_m[m_idx]; dens[n_idx] <- dens_n[n_idx] }
  w <- c(diff(grid), tail(diff(grid),1)); total <- sum(dens*w); if (total > .Machine$double.eps) dens <- dens/total
  list(s=grid, density=dens, varphi=varphi)
}

#' Density mode
#' @export
mode_of_density <- function(s, d){ s[which.max(d)] }

#' One-step recursive volatility + price PI
#' @export
predict_one_step <- function(G, tQML, par, level=0.95){
  R <- c(G[1], diff(G)); dt <- c(tQML[1], diff(tQML)); X <- 0
  for (i in seq_along(R)){
    s2_i <- exp(par$mu + exp(-par$a1*dt[i])*X - par$gamma*par$lambda*par$K*dt[i])
    Zhat <- R[i]/sqrt(s2_i)
    X <- exp(-par$a1*dt[i]) * X + par$theta * Zhat + par$gamma * abs(Zhat) - par$gamma*par$lambda*par$K*dt[i]
  }
  dt_next <- tail(dt, 1); dens <- vol_prediction_density_piecewise(Xn=X, dt_next=dt_next, par=par)
  s2_next <- mode_of_density(dens$s, dens$density); z <- stats::qnorm((1+level)/2)
  PI <- c(tail(G,1) - z*sqrt(s2_next/par$lambda), tail(G,1) + z*sqrt(s2_next/par$lambda))
  list(sigma_next=s2_next, PI=PI, density=dens)
}
