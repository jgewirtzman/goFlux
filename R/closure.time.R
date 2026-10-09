#' Closure duration from elapsed time
#'
#' Duration of a measurement, in seconds, used as \code{t} in
#' \code{\link[goFlux]{MDF}} and in the kappa-max bound of the HM model:
#' the span of the selected window plus one logging interval,
#' \deqn{t = (\max(Etime) - \min(Etime)) + \Delta t,}
#' where \eqn{\Delta t} is the median of the positive differences between
#' consecutive \code{Etime} values.
#'
#' @param Etime numerical vector; elapsed time (seconds) of the observations
#'              retained for one measurement (\code{flag == 1}).
#'
#' @details
#' Each observation stands for one logging interval, so a window of \eqn{n}
#' equally spaced observations covers \eqn{n \cdot \Delta t} seconds. At 1 Hz
#' this equals the earlier convention \code{max(Etime) + 1}; at coarser
#' intervals that convention under-counted the closure by
#' \eqn{\Delta t - 1} seconds (4 s at 5 s logging, 9 s at 10 s). The duration
#' is never the number of observations (\code{nb.obs}); the two coincide only
#' at 1 Hz.
#'
#' The logging interval is the median of the positive time steps, so duplicate
#' timestamps and isolated gaps do not change it. With fewer than two distinct
#' timestamps the interval defaults to 1 s.
#'
#' @return a numerical value (seconds)
#'
#' @seealso \code{\link[goFlux]{MDF}}, \code{\link[goFlux]{goFlux}}
#'
#' @examples
#' closure.time(0:179)              # 180 s at 1 Hz
#' closure.time(seq(0, 175, by = 5)) # 180 s at 5 s logging
#'
#' @export
#'
closure.time <- function(Etime) {
  if (!is.numeric(Etime)) stop("'Etime' must be numeric")
  Etime <- Etime[is.finite(Etime)]
  if (length(Etime) == 0) stop("'Etime' contains no finite values")
  dt <- diff(sort(Etime))
  dt <- dt[dt > 0]
  interval <- if (length(dt) > 0) median(dt) else 1
  diff(range(Etime)) + interval
}
