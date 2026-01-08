#' Consolidate bid/ask quotes at identical timestamps to mid-quotes
#' @export
tba_midquotes <- function(tba, b, a, tgrid = NULL) {
  stopifnot(length(tba) == length(b), length(b) == length(a))
  mid <- (log(a) + log(b)) / 2
  ord <- order(tba)
  tba <- tba[ord]
  mid <- mid[ord]
  runs <- rle(tba)
  idx <- cumsum(runs$lengths)
  start <- c(1, head(idx, -1) + 1)
  out_t <- numeric(length(runs$lengths))
  out_m <- numeric(length(runs$lengths))
  for (i in seq_along(runs$lengths)) {
    rng <- start[i]:idx[i]
    out_t[i] <- tba[rng[1]]
    out_m[i] <- mean(mid[rng])
  }
  if (is.null(tgrid)) data.frame(time = out_t, mid = out_m) else prev_tick(tgrid, out_t, out_m)
}
