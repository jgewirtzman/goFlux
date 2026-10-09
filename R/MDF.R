#' Minimal Detectable Flux (MDF)
#'
#' The minimal detectable flux based on instrument precision and the
#' measurement time: \eqn{MDF = z \cdot p / t \cdot flux.term}, the
#' precision/time form of Christiansen et al. (2015, section 2.5, who define
#' MDF as analytic precision divided by the enclosure time, i.e. z = 1).
#' \code{conf} sets the multiplier z.
#'
#' @param p numerical; precision of the instrument (same units as measured gas; ex. ppm)
#' @param t numerical; measurement time (enclosure time; seconds)
#' @param flux.term numerical; flux term calculated with the function
#'                  \code{\link[goFlux]{flux.term}}
#' @param conf numerical; optional label of the multiplier z applied to the
#'             precision, the standard-normal quantile
#'             \eqn{z = qnorm(1 - (1 - conf)/2)} (1.96 for 0.95, 1.645 for
#'             0.90). It is the Gaussian quantile for a single concentration
#'             difference: z = 1.96 is a benchmark multiplier, not a
#'             calibrated 95 \% test of a fitted flux (see Details). Default
#'             \code{NULL}: z = 1 (unchanged behaviour).
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
#' expressed as one standard deviation, the default MDF (z = 1) is a one-sigma
#' benchmark. An empirical precision can be estimated from the data with
#' \code{\link[goFlux]{empirical.prec}}; a datasheet precision quoted at 1 s
#' can be compared with a record logged at a coarser interval with
#' \code{\link[goFlux]{spec.at.interval}}, if each logged value is an average.
#' The precision given to \code{\link[goFlux]{goFlux}} also sets kappa-max
#' (\code{\link[goFlux]{k.max}}), and so the curvature allowed to the HM model
#' and the model choice in \code{\link[goFlux]{best.flux}}: changing the
#' precision can change the fluxes, not only the MDF.
#'
#' \strong{What the MDF is.} The MDF converts a chosen concentration change
#' (z times the precision) over the closure into a flux. It is a benchmark:
#' it does not depend on the number of observations, and a flux below it
#' cannot be distinguished from zero by this benchmark, which does not make
#' it zero or noise (keep it signed). Regression-based detection (a slope
#' p-value, or a confidence interval for the same linear model) answers a
#' different question: it gains precision with the number of observations and
#' assumes independent residuals (Cowan et al. 2025). Neither covers
#' instrument drift or leaks, which look like a flux (see
#' \code{\link[goFlux]{qc.flags}}).
#'
#' \strong{Attribution.} The precision/time form is that of Christiansen et
#' al. (2015) with z = 1; the "SD x 3 x t99" of that paper is a
#' gas-chromatograph quantification limit from replicate standards (Corley
#' 2003), not a per-closure multiplier. Wassmann et al. (2018) used an
#' empirical, instrument-level precision (the SD of replicate ambient-air
#' analyses) times a multiplier of 3. Parkin, Venterea and Hargreaves (2012)
#' derived detection limits by simulation; their coefficients are specific to
#' their sampling designs.
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
#' Cowan, N., Levy, P., Tigli, M., Toteva, G., & Drewer, J. (2025).
#' Characterisation of analytical uncertainty in chamber soil flux
#' measurements. \emph{European Journal of Soil Science}, 76(2), e70104.
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
