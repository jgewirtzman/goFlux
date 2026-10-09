#' Minimal Detectable Flux (MDF)
#'
#' The minimal detectable flux based on instrument precision and the
#' measurement time: \eqn{MDF = p / t \cdot flux.term} (Christiansen et al.
#' 2015, section 2.5, who define MDF as analytic precision divided by the
#' enclosure time). Optionally scaled to a stated confidence level with
#' \code{conf}.
#'
#' @param p numerical; precision of the instrument (same units as measured gas; ex. ppm)
#' @param t numerical; measurement time (enclosure time; seconds)
#' @param flux.term numerical; flux term calculated with the function
#'                  \code{\link[goFlux]{flux.term}}
#' @param conf numerical; optional two-sided confidence level (e.g. 0.95).
#'             When given, the precision is scaled by the standard-normal
#'             quantile \eqn{z = qnorm(1 - (1 - conf)/2)} (1.96 for 0.95, the
#'             2-sigma criterion; 1.645 for 0.90), i.e.
#'             \eqn{MDF = z \cdot p / t \cdot flux.term}. This generalises the
#'             fixed multiplier of Wassmann et al. (2018), who used k = 3
#'             ("99 \%") on the standard deviation of replicate ambient-air
#'             analyses; for the alpha = 0.05 detection limit see Parkin,
#'             Venterea and Hargreaves (2012). Default \code{NULL} keeps
#'             \eqn{z = 1} (unchanged behaviour).
#'
#' @details
#' \strong{Measurement time \code{t}.} \code{\link[goFlux]{goFlux}} passes
#' the closure duration from \code{\link[goFlux]{closure.time}}: the span of
#' the retained window plus one logging interval, in seconds computed from
#' \code{POSIX.time}. A closure of 180 s logged every 5 s (36 observations)
#' therefore gives the same \eqn{MDF \cdot t / p} as the same closure logged
#' at 1 Hz (180 observations). \code{t} is never the number of observations
#' (\code{nb.obs}). Earlier versions used \code{max(Etime) + 1}, which equals
#' the new value at 1 Hz and under-counts by \eqn{\Delta t - 1} seconds at
#' coarser intervals. The same \code{t} enters the kappa-max bound of the HM
#' model (\code{\link[goFlux]{k.max}}).
#'
#' \strong{Precision \code{p}.} With a datasheet or empirical precision
#' expressed as one standard deviation, the default MDF is a one-sigma limit.
#' Note that the "SD x 3 x t99" of Christiansen et al. (2015) is a
#' gas-chromatograph method quantification limit from replicate standards
#' (Corley 2003), not a per-closure MDF; it should not be read as a
#' \eqn{3 \cdot t} variant of this function. An empirical precision can be
#' estimated from the data with \code{\link[goFlux]{empirical.prec}}; a
#' datasheet precision quoted at 1 s can be rescaled to the logging interval
#' with \code{\link[goFlux]{spec.at.interval}}.
#'
#' @return a numerical value
#'
#' @seealso \code{\link[goFlux]{flux.term}},
#'          \code{\link[goFlux]{empirical.prec}} to estimate \code{p} from the data,
#'          \code{\link[goFlux]{spec.at.interval}}
#'
#' @references Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015).
#' Comparison of
#' \ifelse{html}{\out{CO<sub>2</sub>}}{\eqn{CO[2]}{ASCII}},
#' \ifelse{html}{\out{CH<sub>4</sub>}}{\eqn{CH[4]}{ASCII}} and
#' \ifelse{html}{\out{N<sub>2</sub>O}}{\eqn{N[2]O}{ASCII}} soil-atmosphere
#' exchange measured in static chambers with cavity ring-down spectroscopy and
#' gas chromatography. \emph{Agricultural and Forest Meteorology}, 211-212,
#' 48-57.
#'
#' Corley, J. (2003). Best practices in establishing detection and
#' quantification limits for pesticide residues in foods. In: \emph{Handbook of
#' Residue Analytical Methods for Agrochemicals}. Wiley.
#'
#' Parkin, T. B., Venterea, R. T., & Hargreaves, S. K. (2012). Calculating the
#' detection limits of chamber-based soil greenhouse gas flux measurements.
#' \emph{Journal of Environmental Quality}, 41, 705-715.
#'
#' Wassmann, R., Alberto, M. C., Tirol-Padre, A., Hoang, N. T., Romasanta, R.,
#' Centeno, C. A., & Sander, B. O. (2018). Increasing sensitivity of methane
#' emission measurements in rice through deployment of 'closed chambers' at
#' nighttime. \emph{PLoS ONE}, 13(2), e0191352.
#'
#' @examples
#' ft <- flux.term(V_L = 7.08, P_kPa = 101.3, A_cm2 = 4254, T_C = 26)
#' MDF(p = 1.1, t = 540, flux.term = ft)              # one-sigma benchmark
#' MDF(p = 1.1, t = 540, flux.term = ft, conf = 0.95) # z = 1.96
#'
#' @export
#'
MDF <- function(p, t, flux.term, conf = NULL) {
  z <- if (is.null(conf)) 1 else {
    if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1)
      stop("'conf' must be a single number between 0 and 1")
    stats::qnorm(1 - (1 - conf) / 2)
  }
  z * (p / t) * flux.term
}
