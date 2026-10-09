#' Datasheet precision rescaled to a logging interval
#'
#' Manufacturers quote instrument precision at 1 s (1 Hz). A record logged
#' every \code{dt_s} seconds averages that noise down, so the fair reference
#' for an empirical precision measured on such a record is
#' \eqn{spec_{1s} / \sqrt{dt_s}} (white noise). Compare the
#' \code{\link[goFlux]{empirical.prec}} of a 5 s or 10 s record with this, not
#' with the 1 s figure; otherwise an analyzer can appear to beat its datasheet
#' when it is several times worse. Not valid once \code{dt_s} exceeds the
#' analyzer's Allan minimum (the averaging time beyond which drift dominates).
#'
#' @param spec_1s numerical; datasheet precision at 1 s (same units as the gas).
#' @param dt_s numerical; logging interval in seconds (e.g. the \code{dt_s}
#'             column returned by \code{\link[goFlux]{empirical.prec}}).
#'
#' @returns numerical; same units as \code{spec_1s}, recycled to the longer of
#'          the two arguments.
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{empirical.prec}}, \code{\link[goFlux]{MDF}}
#'
#' @examples
#' spec.at.interval(0.9, c(1, 5, 10))
#' @export
spec.at.interval <- function(spec_1s, dt_s) {
  if(missing(spec_1s)) stop("'spec_1s' is required")
  if(missing(dt_s)) stop("'dt_s' is required")
  if(!is.numeric(spec_1s)) stop("'spec_1s' must be numeric")
  if(!is.numeric(dt_s)) stop("'dt_s' must be numeric")
  if(any(dt_s <= 0, na.rm = TRUE)) stop("'dt_s' must be > 0")
  spec_1s / sqrt(dt_s)
}
