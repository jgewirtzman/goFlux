#' Flux term
#'
#' The "flux term" corrects the flux estimate for atmospheric pressure, air
#' temperature, the effect of water vapor, the volume of the chamber, and
#' the surface area covered by the chamber. This corrections changes the units
#' of the flux estimate from a concentration (ppm or ppb) per time (seconds)
#' into a molarity (µmol or nmol) per area (m2) per time (seconds).
#'
#' @param V_L numerical; total volume of the system (L): chamber or collar,
#'            tubing and the analyzer's internal volume (see Details)
#' @param P_kPa numerical; atmospheric pressure (kPa)
#' @param A_cm2 numerical; area of the soil surface inside the chamber
#'              \ifelse{html}{\out{(cm<sup>2</sup>)}}{\eqn{(cm^2)}{ASCII}}
#' @param T_C numerical; air temperature before chamber closure (Celsius)
#' @param H2O_mol numerical; water vapor concentration in the air before chamber
#'                closure (mol/mol). Default 0 (no correction).
#'
#' @details
#' Flux term units are
#' \ifelse{html}{\out{mol/m<sup>2</sup>}}{\eqn{mol/m^{2}}{ASCII}}
#'
#' Multiplying a slope in ppm/s by the flux term gives
#' \ifelse{html}{\out{µmol m<sup>-2</sup>s<sup>-1</sup>}}{\eqn{\mu mol m^{-2}s^{-1}}{ASCII}},
#' a slope in ppb/s gives
#' \ifelse{html}{\out{nmol m<sup>-2</sup>s<sup>-1</sup>}}{\eqn{nmol m^{-2}s^{-1}}{ASCII}}.
#'
#' \strong{Total volume.} \code{V_L} is the whole closed loop: chamber (or
#' collar plus chamber), tubing and the analyzer's internal volume. Record
#' each component separately so that it can be checked. An x \% error in the
#' volume or the area is an x \% error in the flux and in the
#' \code{\link[goFlux]{MDF}} (Cowan et al. 2025). For the analyzer, use the
#' manufacturer's total sample volume: e.g. 28
#' \ifelse{html}{\out{cm<sup>3</sup>}}{\eqn{cm^3}{ASCII}} for the LI-COR
#' LI-7810 (a value also used for the ABB/LGR GLA131 microportable analyzers,
#' whose datasheet gives about 25
#' \ifelse{html}{\out{cm<sup>3</sup>}}{\eqn{cm^3}{ASCII}}). The 70
#' \ifelse{html}{\out{cm<sup>3</sup>}}{\eqn{cm^3}{ASCII}} in the
#' instrument table of the goFlux website (\code{example_auxfile.xlsx}) is
#' for the larger UGGA, not for the microportable GLA131 ("MGGA").
#'
#' @return a numerical value
#'
#' @references
#' Cowan, N., Levy, P., Tigli, M., Toteva, G., & Drewer, J. (2025).
#' Characterisation of analytical uncertainty in chamber soil flux
#' measurements. \emph{European Journal of Soil Science}, 76(2), e70104.
#' \doi{10.1111/ejss.70104}
#'
#' @examples
#' # 1 L chamber on a 100 cm2 footprint, 20 C, 101.3 kPa
#' ft <- flux.term(1, 101.3, 100, 20)
#' 0.05 * ft   # a CH4 slope of 0.05 ppb/s, in nmol m-2 s-1
#'
#' @export
#'
flux.term <- function(V_L, P_kPa, A_cm2, T_C, H2O_mol = 0) {
  (V_L * P_kPa * (1 - H2O_mol)) / (8.314 * (A_cm2/10000) * (T_C + 273.15))
}
