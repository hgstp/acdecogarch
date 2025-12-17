
#' Core QML estimator for ECOGARCH(1,1) on irregularly spaced data
#' @param Gdt Returns at jump times
#' @param dt Interarrival times
#' @param init list(a1, theta, gamma, mu)
#' @param method recursion type ('approx'|'exact')
#' @param grad gradient option ('none'|'fd')
#' @param h FD step size
#' @export
ecogarch_qmle <- function(Gdt, dt, init=list(a1=0.1, theta=-0.1, gamma=0.3, mu=-3),
                             method=c('approx','exact'), grad=c('none','fd'), h=1e-4){
  method <- match.arg(method); grad <- match.arg(grad);
  stopifnot(length(Gdt) == length(dt))
  n <- length(Gdt)
  lam_hat <- n / sum(dt)
  K_hat <- sqrt(2/(pi*lam_hat))
  forward <- function(par){
    a1 <- exp(par[1]); theta <- par[2]; gamma <- par[3]; mu <- par[4]
    X <- numeric(n); s2 <- numeric(n); Zhat <- numeric(n)
    for (i in seq_len(n)){
      Xprev <- if (i==1) 0 else X[i-1]
      if (method == 'approx'){
        s2[i] <- exp(mu + exp(-a1*dt[i]) * Xprev - gamma*lam_hat*K_hat*dt[i])
        Zhat[i] <- Gdt[i] / sqrt(s2[i])
        X[i] <- exp(-a1*dt[i]) * Xprev + theta * Zhat[i] + gamma * abs(Zhat[i]) - gamma*lam_hat*K_hat*dt[i]
      } else {
        s2[i] <- exp(mu + exp(-a1*dt[i]) * Xprev - gamma*lam_hat*K_hat*(1-exp(-a1*dt[i]))/a1)
        Zhat[i] <- Gdt[i] / sqrt(s2[i])
        X[i] <- exp(-a1*dt[i]) * Xprev + theta*Zhat[i] + gamma*(abs(Zhat[i]) - lam_hat*K_hat*(1-exp(-a1*dt[i]))/a1)
      }
    }
    list(a1=a1, theta=theta, gamma=gamma, mu=mu, X=X, s2=s2, Zhat=Zhat)
  }
  qll <- function(par){
      ff <- forward(par)
      v <- ff$s2/lam_hat
      0.5*sum(log(v) + (Gdt^2)/v)
  }

  gr_fd <- function(par){
    p <- length(par)
    g <- numeric(p)
    for(k in seq_len(p)){
      step <- rep(0,p)
      step[k] <- h
      g[k] <- (qll(par+step)-qll(par-step))/(2*h)
    }
    g
  }

  par0 <- c(log(init$a1), init$theta, init$gamma, init$mu)
  lower <- c(log(1e-6), -Inf, -Inf, -Inf); upper <- c(log(10), Inf, Inf, Inf)
  fit <- if (grad == 'fd') stats::optim(par0, qll, gr=gr_fd, method='L-BFGS-B', lower=lower, upper=upper)
         else               stats::optim(par0, qll,              method='L-BFGS-B', lower=lower, upper=upper)
  ff <- forward(fit$par)
  list(par=list(a1=ff$a1, theta=ff$theta, gamma=ff$gamma, mu=ff$mu, lambda=lam_hat, K=K_hat),
       X=ff$X, s2=ff$s2, Zhat=ff$Zhat,
       value=fit$value, convergence=fit$convergence, counts=fit$counts, method=method, grad=grad, h=h)
}
