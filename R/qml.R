
#' QML wrappers (mu and HF variants)
#' @export
ecogarch_qmle_mu <- function(G, tQML, init=list(a1=0.1, theta=-0.1, gamma=0.3, mu=-3)){
  R <- c(G[1], diff(G)); dt <- c(tQML[1], diff(tQML)); ecogarch_qmle(R, dt, init=init, method='approx') }
#' @export
ecogarch_qmle_hfdata <- function(G, tQML, init=list(a1=0.1, theta=-0.1, gamma=0.3, mu=-3), lag_offset=10){
  if (length(G) <= 2*lag_offset+1) stop('Not enough observations for HF variant')
  R <- G[(lag_offset+1):(length(G)-lag_offset)] - G[lag_offset:(length(G)-lag_offset-1)]
  dt <- tQML[(lag_offset+1):(length(tQML)-lag_offset)] - tQML[lag_offset:(length(tQML)-lag_offset-1)]
  ecogarch_qmle(R, dt, init=init, method='approx') }
