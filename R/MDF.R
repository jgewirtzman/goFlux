#' Minimal Detectable Flux (MDF)
#'
#' The minimal detectable flux based on instrument precision and the
#' measurement time.
#'
#' @param p numerical; precision of the instrument (same units as measured gas; ex. ppm)
#' @param t numerical; measurement time (enclosure time; seconds)
#' @param flux.term numerical; flux term calculated with the function
#'                  \code{\link[goFlux]{flux.term}}
#' @param conf numerical; optional two-sided confidence level (e.g. 0.95).
#'             When given, the precision is scaled by the standard-normal
#'             quantile \eqn{z = qnorm(1 - (1 - conf)/2)} (1.96 for 0.95, the
#'             2-sigma criterion; 1.645 for 0.90), i.e.
#'             \eqn{MDF = z \cdot p / t \cdot flux.term} (Wassmann et al. 2018).
#'             Default \code{NULL} keeps \eqn{z = 1} (unchanged behaviour).
#'
#' @return a numerical value
#'
#' @seealso \code{\link[goFlux]{flux.term}}
#'
#' @references Christiansen et al. (2015). Comparison of
#' \ifelse{html}{\out{CO<sub>2</sub>}}{\eqn{CO[2]}{ASCII}},
#' \ifelse{html}{\out{CH<sub>4</sub>}}{\eqn{CH[4]}{ASCII}} and
#' \ifelse{html}{\out{N<sub>2</sub>O}}{\eqn{N[2]O}{ASCII}} soil-atmosphere
#' exchange measured in static chambers with cavity ring-down spectroscopy and
#' gas chromatography. \emph{Agricultural and Forest Meteorology}, 211, 48-57.
#'
#' Wassmann, R., et al. (2018). Increasing sensitivity of methane emission
#' measurements in rice through deployment of 'closed chambers' at nighttime.
#' \emph{PLoS ONE}, 13(2), e0191352.
#'
#' @seealso \code{\link[goFlux]{empirical.prec}} to estimate \code{p} from the data.
#'
#' @keywords internal
#'
MDF <- function(p, t, flux.term, conf = NULL) {
  z <- if (is.null(conf)) 1 else {
    if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1)
      stop("'conf' must be a single number between 0 and 1")
    stats::qnorm(1 - (1 - conf) / 2)
  }
  z * (p / t) * flux.term
}
