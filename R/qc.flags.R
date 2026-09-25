#' Post-hoc quality flags for flux results
#'
#' Adds simple, independent logical flags to the output of
#' \code{\link[goFlux]{best.flux}} (or \code{\link[goFlux]{goFlux}}). Nothing is
#' removed or modified: each flag is a new column (\code{TRUE} = the check
#' fired, \code{NA} = could not be evaluated), and \code{qc.any} says whether
#' any fired. Filter on them downstream as you see fit.
#'
#' @param flux.result data.frame; output from \code{\link[goFlux]{best.flux}}
#'                    or \code{\link[goFlux]{goFlux}}.
#' @param dataframe data.frame; optional, the flagged concentration data the
#'                  fluxes were computed from (output of
#'                  \code{\link[goFlux]{click.peak2}} or
#'                  \code{\link[goFlux]{crop.meas}}). Needed for
#'                  \code{qc.ambient} and \code{qc.noisy}.
#' @param gastype character string; the gas column of \code{dataframe}.
#' @param c0.mult numerical; \code{qc.c0} fires when \code{C0} exceeds
#'                \code{c0.mult} times the median \code{C0} of its group
#'                (\code{by}): the chamber was not at ambient at closure.
#'                Default 1.5. \code{NULL} skips the flag.
#' @param by character string; optional column of \code{flux.result} defining
#'           the groups (campaign, day, instrument) for \code{qc.c0} and
#'           \code{qc.noisy}. Default \code{NULL}: all rows.
#' @param min.obs numerical; \code{qc.min.obs} fires when \code{nb.obs} is below
#'                it. Default 60. \code{NULL} skips the flag.
#' @param ambient.sigma,ambient.secs numerical; \code{qc.ambient} fires when
#'                the mean of the first \code{ambient.secs} seconds of the
#'                flagged window differs from the pre-closure ambient (rows of
#'                \code{dataframe} before \code{start.time}) by more than
#'                \code{ambient.sigma} times the closure's Allan precision:
#'                the chamber was already enriched when the window started.
#'                Defaults 3 and 10. Needs \code{dataframe}.
#' @param noisy.mult numerical; \code{qc.noisy} fires when the closure's own
#'                precision (MAD of first differences / sqrt(2) on its flagged
#'                rows) exceeds \code{noisy.mult} times the record (or group)
#'                precision from \code{\link[goFlux]{empirical.prec}}: a
#'                disturbed closure. When the record was logged at more than
#'                one interval, each closure is compared with the record
#'                precision at its own logging interval (the \code{dt_s} of
#'                \code{empirical.prec}). Default 1.5. Needs \code{dataframe}.
#'
#' @details
#' \code{qc.convex} fires when the Hutchinson-Mosier curvature is convex
#' (\code{HM.k < 0}, accelerating concentration change, which no chamber
#' mechanism produces and which points at a leak, a disturbance or a delayed
#' seal); it can only be \code{TRUE} when \code{\link[goFlux]{goFlux}} was run
#' with \code{k.min < 0}, otherwise it is \code{FALSE} for a fitted HM and
#' \code{NA} when \code{HM.k} is missing.
#'
#' The CO2 tracer test (CO2 must accumulate in a sealed chamber on a respiring
#' surface) is not a flag here because it is a plain comparison between two
#' \code{best.flux} outputs: see \code{\link[goFlux]{co2.tracer}}.
#'
#' @returns \code{flux.result} with the columns \code{qc.c0}, \code{qc.convex},
#'          \code{qc.min.obs}, \code{qc.ambient}, \code{qc.noisy} (those
#'          requested) and \code{qc.any} appended, plus the helper values
#'          \code{qc.c0.ratio}, \code{qc.ambient.dev} (in sigma units) and
#'          \code{qc.noisy.ratio}.
#'
#' @include goFlux-package.R
#' @include empirical.prec.R
#'
#' @seealso \code{\link[goFlux]{co2.tracer}}, \code{\link[goFlux]{best.flux}}
#'
#' @examples
#' data(manID.UGGA)
#' CH4_flux <- goFlux(manID.UGGA, "CH4dry_ppb")
#' CH4_best <- best.flux(CH4_flux)
#' CH4_qc <- qc.flags(CH4_best, dataframe = manID.UGGA, gastype = "CH4dry_ppb")
#' CH4_qc[, c("UniqueID", "best.flux", "qc.c0", "qc.convex", "qc.min.obs",
#'            "qc.ambient", "qc.noisy", "qc.any")]
#' @export
qc.flags <- function(flux.result, dataframe = NULL, gastype = NULL,
                     c0.mult = 1.5, by = NULL, min.obs = 60,
                     ambient.sigma = 3, ambient.secs = 10, noisy.mult = 1.5) {

  # Check arguments
  if(missing(flux.result)) stop("'flux.result' is required")
  if(!is.data.frame(flux.result)) stop("'flux.result' must be of class data.frame")
  if(!any(grepl("\\<UniqueID\\>", names(flux.result)))) stop("'flux.result' must contain 'UniqueID'")
  if(!is.null(dataframe)){
    if(!is.data.frame(dataframe)) stop("'dataframe' must be of class data.frame")
    if(is.null(gastype)) stop("'gastype' is required when 'dataframe' is provided")
    for(col in c("UniqueID", "flag", "Etime", gastype)){
      if(!any(grepl(paste("\\<", col, "\\>", sep = ""), names(dataframe)))){
        stop("'dataframe' must contain the column '", col, "'")}}}
  if(!is.null(by) && !any(grepl(paste("\\<", by, "\\>", sep = ""), names(flux.result)))){
    stop("'flux.result' must contain a column that matches 'by'")}

  fx <- as.data.frame(flux.result)
  n <- nrow(fx)
  uid <- as.character(fx$UniqueID)
  g <- if(is.null(by)) rep("all", n) else as.character(fx[[by]])
  flags <- list()

  # qc.c0: starting concentration vs group median
  if(!is.null(c0.mult)){
    if(any(grepl("\\<C0\\>", names(fx)))){
      med <- tapply(fx$C0, g, stats::median, na.rm = TRUE)
      fx$qc.c0.ratio <- as.numeric(fx$C0) / as.numeric(med[g])
      fx$qc.c0 <- ifelse(is.finite(fx$qc.c0.ratio), fx$qc.c0.ratio > c0.mult, NA)
    } else fx$qc.c0 <- NA
    flags <- c(flags, "qc.c0")
  }

  # qc.convex: accelerating HM curvature
  fx$qc.convex <- if(any(grepl("\\<HM.k\\>", names(fx)))) ifelse(is.na(fx$HM.k), NA, fx$HM.k < 0) else NA
  flags <- c(flags, "qc.convex")

  # qc.min.obs
  if(!is.null(min.obs)){
    fx$qc.min.obs <- if(any(grepl("\\<nb.obs\\>", names(fx)))) ifelse(is.na(fx$nb.obs), NA, fx$nb.obs < min.obs) else NA
    flags <- c(flags, "qc.min.obs")
  }

  # flags needing the concentration data
  if(!is.null(dataframe)){
    d <- dataframe[!is.na(dataframe$flag) & dataframe$flag == 1, ]
    sp <- split(d, as.character(d$UniqueID))
    per <- function(f){ r <- vapply(sp, f, numeric(1)); r[match(uid, names(sp))] }

    # qc.ambient: window start vs pre-closure ambient
    if(!is.null(ambient.sigma) && any(grepl("\\<start.time\\>", names(dataframe)))){
      pre <- dataframe[dataframe$POSIX.time < dataframe$start.time, ]
      amb <- tapply(pre[[gastype]], as.character(pre$UniqueID), mean, na.rm = TRUE)
      dev <- per(function(x){
        a <- amb[as.character(x$UniqueID[1])]
        if(is.null(a) || is.na(a)) return(NA_real_)
        v <- x[[gastype]][!is.na(x[[gastype]])]
        s <- stats::sd(diff(v)) / sqrt(2)
        abs(mean(x[[gastype]][x$Etime <= ambient.secs], na.rm = TRUE) - a) / s })
      fx$qc.ambient.dev <- dev
      fx$qc.ambient <- ifelse(is.na(dev), NA, dev > ambient.sigma)
      flags <- c(flags, "qc.ambient")
    }

    # qc.noisy: closure precision vs record (group) precision
    if(!is.null(noisy.mult)){
      s.clo <- per(function(x){ v <- x[[gastype]][!is.na(x[[gastype]])]
        if(length(v) < 3) NA_real_ else stats::mad(diff(v), constant = 1.4826) / sqrt(2) })
      # empirical.prec() returns one row per logging interval (dt_s) when the
      # record changed interval. Compare each closure with the record (group)
      # precision at the closure's own interval (median spacing of its flagged
      # rows, nearest 0.5 s); fall back to the interval with the most
      # differences when none matches (or without POSIX.time).
      dt.clo <- if(any(grepl("\\<POSIX.time\\>", names(dataframe)))){
        per(function(x){ dd <- diff(as.numeric(x$POSIX.time)); dd <- dd[is.finite(dd) & dd > 0]
          if(length(dd) < 1) NA_real_ else round(stats::median(dd) * 2) / 2 })
      } else rep(NA_real_, n)
      pick <- function(ep, dt){
        if(nrow(ep) == 0) return(NA_real_)
        i <- if(is.na(dt)) integer(0) else which(!is.na(ep$dt_s) & ep$dt_s == dt)
        if(length(i) == 0) i <- which.max(ifelse(is.na(ep$n), -1, ep$n))
        ep$prec[i[1]] }
      if(is.null(by)){
        ep <- empirical.prec(dataframe, gastype)
        s.rec <- vapply(dt.clo, function(dt) pick(ep, dt), numeric(1))
      } else {
        tg <- g[match(as.character(dataframe$UniqueID), uid)]
        dd <- dataframe; dd$qc_by <- tg
        ep <- empirical.prec(dd[!is.na(dd$qc_by), ], gastype, by = "qc_by")
        s.rec <- vapply(seq_len(n), function(i)
          pick(ep[ep$qc_by == g[i], , drop = FALSE], dt.clo[i]), numeric(1))
      }
      fx$qc.noisy.ratio <- s.clo / s.rec
      fx$qc.noisy <- ifelse(is.finite(fx$qc.noisy.ratio), fx$qc.noisy.ratio > noisy.mult, NA)
      flags <- c(flags, "qc.noisy")
    }
  }

  m <- as.matrix(fx[, unlist(flags), drop = FALSE])
  fx$qc.any <- apply(m, 1, function(r) if(all(is.na(r))) NA else any(r, na.rm = TRUE))
  fx
}

#' CO2 tracer check between two flux results
#'
#' On a respiring surface a sealed chamber must accumulate CO2: a CO2 flux
#' that is not clearly positive means the chamber leaked, was not sealed, or
#' was on a surface without respiratory CO2 efflux. This is a plain
#' post-hoc comparison of two \code{\link[goFlux]{best.flux}} outputs (CO2 and
#' the gas of interest), joined on \code{UniqueID}.
#'
#' @param co2.flux.result data.frame; \code{\link[goFlux]{best.flux}} output for
#'                        \code{"CO2dry_ppm"}.
#' @param flux.result data.frame; \code{\link[goFlux]{best.flux}} output for the
#'                    gas of interest (e.g. CH4), whose row order is returned.
#' @param p.val numerical; significance level for the linear-model p-value of
#'              the CO2 flux. Default 0.05.
#'
#' @returns A logical vector aligned with the rows of \code{flux.result}:
#'          \code{TRUE} when the CO2 \code{best.flux} is positive and its
#'          \code{LM.p.val} is below \code{p.val}; \code{NA} when the
#'          \code{UniqueID} has no CO2 result. Do not use it on open water,
#'          lit foliage or dead wood.
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{qc.flags}}
#'
#' @examples
#' data(manID.UGGA)
#' CO2_best <- best.flux(goFlux(manID.UGGA, "CO2dry_ppm"))
#' CH4_best <- best.flux(goFlux(manID.UGGA, "CH4dry_ppb"))
#' CH4_best$co2.tracer <- co2.tracer(CO2_best, CH4_best)
#' CH4_best[, c("UniqueID", "best.flux", "co2.tracer")]
#' @export
co2.tracer <- function(co2.flux.result, flux.result, p.val = 0.05) {
  if(!is.data.frame(co2.flux.result) | !is.data.frame(flux.result))
    stop("'co2.flux.result' and 'flux.result' must be of class data.frame")
  for(col in c("UniqueID", "best.flux", "LM.p.val")){
    if(!any(grepl(paste("\\<", col, "\\>", sep = ""), names(co2.flux.result)))){
      stop("'co2.flux.result' must contain the column '", col, "'")}}
  i <- match(as.character(flux.result$UniqueID), as.character(co2.flux.result$UniqueID))
  ok <- co2.flux.result$best.flux > 0 & co2.flux.result$LM.p.val < p.val
  ifelse(is.na(i), NA, ok[i])
}
