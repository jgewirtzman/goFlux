#' Detection class of each flux relative to the minimal detectable flux
#'
#' Classifies each flux of a \code{\link[goFlux]{best.flux}} (or
#' \code{\link[goFlux]{goFlux}}) table as \code{"emission"}, \code{"uptake"} or
#' \code{"below MDF"} by comparing it with the minimal detectable flux
#' (\code{\link[goFlux]{MDF}}). Nothing is removed or modified: the precision,
#' duration, MDF and class used are appended as new columns, and fluxes below
#' the MDF keep their sign and value.
#'
#' @param flux.result data.frame; output from \code{\link[goFlux]{best.flux}}
#'                    or \code{\link[goFlux]{goFlux}}. Must contain
#'                    \code{UniqueID}, \code{flux.term} and \code{flux.col}.
#' @param dataframe data.frame; optional, the flagged concentration data the
#'                  fluxes were computed from (output of
#'                  \code{\link[goFlux]{click.peak2}} or
#'                  \code{\link[goFlux]{crop.meas}}). Used for the closure
#'                  duration (\code{\link[goFlux]{closure.time}} of the rows
#'                  with \code{flag == 1}) and, when \code{prec = NULL}, for the
#'                  empirical precision.
#' @param gastype character string; the gas column of \code{dataframe}.
#'                Required with \code{dataframe}.
#' @param prec precision (1 sigma, units of \code{gastype}) on which the MDF is
#'             based. One of: \code{NULL} (default), see Details; a single
#'             number; a numeric vector aligned with the rows of
#'             \code{flux.result}; or a data.frame with the columns
#'             \code{UniqueID} and \code{prec} (e.g. per-closure values, or a
#'             group value repeated for each of its closures).
#' @param by character string; optional column of \code{flux.result} defining
#'           the groups (e.g. analyzer x field day) over which the empirical
#'           precision is summarised when \code{prec = NULL}. Default
#'           \code{NULL}: all rows form one group.
#' @param conf numerical; label of the multiplier applied to the precision,
#'             \eqn{z = qnorm(1 - (1 - conf)/2)} (see \code{\link[goFlux]{MDF}}).
#'             Default 0.95 (z = 1.96). \code{NULL} gives z = 1.
#' @param t numerical; optional closure duration in seconds, a single number
#'          or a vector aligned with the rows of \code{flux.result}. Default
#'          \code{NULL}: from \code{dataframe} with
#'          \code{\link[goFlux]{closure.time}}.
#' @param flux.col character string; the flux column to classify. Default
#'                 \code{"best.flux"}.
#'
#' @details
#' The MDF of each row is \eqn{MDF = z \cdot prec / t \cdot flux.term}
#' (\code{\link[goFlux]{MDF}}), with \code{t} the closure duration in seconds
#' (span of the window plus one logging interval), never \code{nb.obs}. The
#' precision is, in order:
#' \itemize{
#'   \item \code{prec}, when given (a datasheet value, a bench sigma, a sigma
#'         computed elsewhere);
#'   \item with \code{prec = NULL} and \code{dataframe}: the group precision of
#'         \code{\link[goFlux]{empirical.prec}} (\code{method = "hadamard"}),
#'         the median over the closures of each \code{by} group of the
#'         per-closure second-difference precision, at the closure's own
#'         logging interval;
#'   \item with \code{prec = NULL} and no \code{dataframe}: goFlux's own
#'         \code{MDF} column is used as it is (it was computed with the
#'         \code{prec} and \code{conf} given to \code{\link[goFlux]{goFlux}};
#'         \code{conf} and \code{t} here are then ignored).
#' }
#' A flux is \code{"emission"} when \code{flux > MDF}, \code{"uptake"} when
#' \code{flux < -MDF}, and \code{"below MDF"} otherwise.
#'
#' \strong{What the class means.} The MDF is a benchmark built on a single
#' concentration difference: \code{conf = 0.95} labels the Gaussian quantile
#' 1.96 applied to the precision, not a calibrated 95 \% test of the fitted
#' flux. "Below MDF" means that the flux cannot be distinguished from zero by
#' this benchmark; it does not mean that the flux is zero or noise, and the
#' value should stay signed in means and budgets. Regression-based detection
#' (slope p-value, or a confidence interval for the same linear model) answers
#' a different question: it gains precision with the number of observations
#' and assumes independent residuals (Cowan et al. 2025). A linear-model
#' p-value is not the test of a Hutchinson-Mosier estimate.
#'
#' The precision given to \code{\link[goFlux]{goFlux}} (\code{prec} or the
#' \code{*_prec} columns) also sets kappa-max and so the curvature allowed to
#' the HM model and the model choice of \code{\link[goFlux]{best.flux}}.
#' \code{flux.class()} only reclassifies the fluxes it is given; to see the
#' effect of a precision on the fluxes themselves, refit with that precision.
#'
#' Neither the precision nor the MDF covers instrument drift or leaks, which
#' look like a flux; see \code{\link[goFlux]{qc.flags}} (\code{leak.rate},
#' \code{blank.slope}).
#'
#' @returns \code{flux.result} with the columns \code{det.prec} (precision
#'          used), \code{det.t} (duration, s), \code{det.MDF} (MDF, flux
#'          units) and \code{det.class} appended. With \code{prec = NULL} and
#'          \code{dataframe}, also \code{det.prec.closure}, the closure's own
#'          second-difference precision (a diagnostic, not the MDF basis).
#'
#' @references
#' Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015). Comparison of
#' \ifelse{html}{\out{CO<sub>2</sub>}}{\eqn{CO[2]}{ASCII}},
#' \ifelse{html}{\out{CH<sub>4</sub>}}{\eqn{CH[4]}{ASCII}} and
#' \ifelse{html}{\out{N<sub>2</sub>O}}{\eqn{N[2]O}{ASCII}} soil-atmosphere
#' exchange measured in static chambers with cavity ring-down spectroscopy and
#' gas chromatography. \emph{Agricultural and Forest Meteorology}, 211-212,
#' 48-57. \doi{10.1016/j.agrformet.2015.06.004}
#'
#' Cowan, N., Levy, P., Tigli, M., Toteva, G., & Drewer, J. (2025).
#' Characterisation of analytical uncertainty in chamber soil flux
#' measurements. \emph{European Journal of Soil Science}, 76(2), e70104.
#' \doi{10.1111/ejss.70104}
#'
#' Wassmann, R., Alberto, M. C., Tirol-Padre, A., Hoang, N. T., Romasanta, R.,
#' Centeno, C. A., & Sander, B. O. (2018). Increasing sensitivity of methane
#' emission measurements in rice through deployment of 'closed chambers' at
#' nighttime. \emph{PLoS ONE}, 13(2), e0191352.
#' \doi{10.1371/journal.pone.0191352}
#'
#' @include goFlux-package.R
#' @include empirical.prec.R
#'
#' @seealso \code{\link[goFlux]{MDF}}, \code{\link[goFlux]{empirical.prec}},
#'          \code{\link[goFlux]{closure.time}}, \code{\link[goFlux]{qc.flags}}
#'
#' @examples
#' data(manID.UGGA)
#' CH4_best <- best.flux(goFlux(manID.UGGA, "CH4dry_ppb"))
#' # empirical (second-difference) precision of the closures, z = 1.96
#' CH4_det <- flux.class(CH4_best, dataframe = manID.UGGA, gastype = "CH4dry_ppb")
#' CH4_det[, c("UniqueID", "best.flux", "det.prec", "det.t", "det.MDF", "det.class")]
#' # a precision computed elsewhere
#' flux.class(CH4_best, dataframe = manID.UGGA, gastype = "CH4dry_ppb", prec = 1.2)$det.MDF
#' @export
flux.class <- function(flux.result, dataframe = NULL, gastype = NULL, prec = NULL,
                       by = NULL, conf = 0.95, t = NULL, flux.col = "best.flux") {

  # Check arguments
  if(missing(flux.result)) stop("'flux.result' is required")
  if(!is.data.frame(flux.result)) stop("'flux.result' must be of class data.frame")
  for(col in c("UniqueID", "flux.term", flux.col)){
    if(!col %in% names(flux.result)) stop("'flux.result' must contain the column '", col, "'")}
  if(!is.null(dataframe)){
    if(!is.data.frame(dataframe)) stop("'dataframe' must be of class data.frame")
    if(is.null(gastype)) stop("'gastype' is required when 'dataframe' is provided")
    for(col in c("UniqueID", "flag", gastype)){
      if(!col %in% names(dataframe)) stop("'dataframe' must contain the column '", col, "'")}
    if(!any(c("Etime", "POSIX.time") %in% names(dataframe))){
      stop("'dataframe' must contain 'Etime' or 'POSIX.time'")}}
  if(!is.null(by) && !by %in% names(flux.result)){
    stop("'flux.result' must contain a column that matches 'by'")}

  fx <- as.data.frame(flux.result)
  n <- nrow(fx)
  uid <- as.character(fx$UniqueID)
  flux <- fx[[flux.col]]

  if(is.null(prec) && is.null(dataframe)){
    # goFlux's own MDF, as computed with the prec and conf given to goFlux()
    if(!"MDF" %in% names(fx)) stop("'flux.result' has no 'MDF' column; give 'prec' or 'dataframe'")
    fx$det.prec <- if("prec" %in% names(fx)) fx$prec else NA_real_
    fx$det.t <- NA_real_
    m <- abs(fx$MDF)
  } else {
    # duration
    if(is.null(t)){
      if(is.null(dataframe)) stop("give 'dataframe' (for the closure duration) or 't'")
      dd <- as.data.frame(dataframe)
      dd <- dd[!is.na(dd$flag) & dd$flag == 1, ]
      et <- if("POSIX.time" %in% names(dd)) as.numeric(dd$POSIX.time) else dd$Etime
      tt <- vapply(split(et, as.character(dd$UniqueID)), closure.time, numeric(1))
      t <- unname(tt[uid])
    } else {
      if(!is.numeric(t) || !(length(t) %in% c(1, n))) stop("'t' must be numeric, of length 1 or nrow(flux.result)")
      t <- rep_len(t, n)
    }
    # precision
    if(is.null(prec)){
      g <- if(is.null(by)) rep("all", n) else as.character(fx[[by]])
      pc <- prec.by.closure(dataframe, gastype, uid, g)
      p <- pc$group
      fx$det.prec.closure <- pc$closure
    } else if(is.data.frame(prec)){
      if(!all(c("UniqueID", "prec") %in% names(prec))) stop("a data.frame 'prec' must contain 'UniqueID' and 'prec'")
      p <- prec$prec[match(uid, as.character(prec$UniqueID))]
    } else {
      if(!is.numeric(prec) || !(length(prec) %in% c(1, n))) stop("'prec' must be numeric, of length 1 or nrow(flux.result), or a data.frame")
      p <- rep_len(prec, n)
    }
    fx$det.prec <- p
    fx$det.t <- t
    m <- abs(MDF(p, t, fx$flux.term, conf = conf))
  }
  fx$det.MDF <- m
  fx$det.class <- ifelse(is.na(m) | is.na(flux), NA_character_,
                         ifelse(flux > m, "emission", ifelse(flux < -m, "uptake", "below MDF")))
  fx
}
