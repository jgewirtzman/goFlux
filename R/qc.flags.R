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
#'                  \code{\link[goFlux]{crop.meas}}), including the rows
#'                  before the closure start (the \code{\link[goFlux]{obs.win}}
#'                  shoulder). Needed for \code{qc.ambient}, \code{qc.clock}
#'                  and \code{qc.noisy}.
#' @param gastype character string; the gas column of \code{dataframe}.
#' @param c0.mult numerical; \code{qc.c0} fires when \code{C0} exceeds
#'                \code{c0.mult} times the median \code{C0} of its group
#'                (\code{by}): the chamber was not at ambient at closure.
#'                Default 1.5. \code{NULL} skips the flag.
#' @param by character string; optional column of \code{flux.result} defining
#'           the groups (campaign, day, instrument) for \code{qc.c0} and for the
#'           group precision used by \code{qc.ambient} and \code{qc.noisy}.
#'           Default \code{NULL}: all rows.
#' @param min.obs numerical; \code{qc.min.obs} fires when \code{nb.obs} is below
#'                it. Default 60. \code{NULL} skips the flag.
#' @param ambient.sigma numerical; \code{qc.ambient} fires when the headspace
#'                at the moment of sealing differs from the pre-closure ambient
#'                by more than \code{ambient.sigma} times the tolerance (see
#'                Details). Default 3. \code{NULL} skips \code{qc.ambient} and
#'                \code{qc.clock}. Needs \code{dataframe}.
#' @param ambient.secs numerical; length (s) of the sealing interval after the
#'                recorded closure start whose median is compared with ambient.
#'                Default 10.
#' @param ambient.pre numerical; length (s) of the pre-closure record before
#'                the recorded closure start that gives the ambient median and
#'                its spread. Default 60.
#' @param seal.time character string; the column of \code{dataframe} holding
#'                the recorded closure start (the moment of sealing). Default
#'                \code{"start.time"}, the auxfile start carried by
#'                \code{\link[goFlux]{obs.win}} and
#'                \code{\link[goFlux]{click.peak2}}. Do not use
#'                \code{start.time_corr} (the start of the fitting window) or
#'                the \code{start.time} returned by
#'                \code{\link[goFlux]{crop.meas}}, which has been moved by the
#'                dead band: keep the original start in a column of its own and
#'                name it here.
#' @param noisy.mult numerical; \code{qc.noisy} fires when the closure's own
#'                precision exceeds \code{noisy.mult} times the group precision
#'                (both from \code{\link[goFlux]{empirical.prec}}, second
#'                differences of the flagged rows; the group value is the median
#'                over its closures, at the closure's own logging interval): a
#'                disturbed closure. Default 1.5. Needs \code{dataframe}.
#'
#' @details
#' \strong{Ambient start (\code{qc.ambient}).} The check is made at the moment
#' of sealing, not at the start of the fitting window, which follows the
#' dead-band transient and may legitimately be enriched already. The median of
#' the first \code{ambient.secs} seconds after the recorded closure start
#' (\code{seal.time}) is compared with the median of the \code{ambient.pre}
#' seconds before it. The flag fires when
#' \deqn{|level_{seal} - ambient| > ambient.sigma \times \max(\sigma, MAD_{pre})}{|seal - ambient| > ambient.sigma x max(sigma, MAD_pre)}
#' where \eqn{\sigma} is the group precision from
#' \code{\link[goFlux]{empirical.prec}} and \eqn{MAD_{pre}}{MAD_pre} the
#' normal-consistent MAD of the pre-closure record (ambient air varies more
#' than the analyzer noise). \code{qc.ambient.dev} is the signed difference in
#' units of that tolerance. The pre-closure trend over \code{ambient.pre}
#' (\code{qc.ambient.drift}, in the same units) is reported but never raises
#' the flag: air outside a sealed chamber cannot bias its slope. A closure that
#' rises steeply moves measurably within \code{ambient.secs}; shorten it for
#' fast-rising closures.
#'
#' \strong{Clock mismatch (\code{qc.clock}).} When the fitting window begins
#' more than one logging interval before the recorded start, the field and
#' analyzer clocks disagree and the "pre-closure" record may lie inside the
#' closure: \code{qc.clock} is \code{TRUE} and \code{qc.ambient} is \code{NA}
#' (not evaluated).
#'
#' \strong{Drift is not flagged here.} A difference-based precision is blind
#' to instrument drift and to slow leaks, which change the concentration
#' smoothly over a closure and look like a flux. None of these flags measures
#' them. Make periodic blank closures on an inert surface (a sealed chamber on
#' a plate or foil): their apparent flux is the drift or leak floor to compare
#' fluxes with. The ambient-shoulder trend (\code{qc.ambient.drift}) can flag
#' an unstable analyzer but cannot measure drift inside the chamber.
#'
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
#'          \code{qc.min.obs}, \code{qc.ambient}, \code{qc.clock},
#'          \code{qc.noisy} (those requested) and \code{qc.any} appended, plus
#'          the helper values \code{qc.c0.ratio}, \code{qc.ambient.dev} and
#'          \code{qc.ambient.drift} (signed, in units of the ambient
#'          tolerance), \code{qc.prec} (the group precision used) and
#'          \code{qc.noisy.ratio}.
#'
#' @include goFlux-package.R
#' @include empirical.prec.R
#'
#' @seealso \code{\link[goFlux]{co2.tracer}}, \code{\link[goFlux]{best.flux}},
#'          \code{\link[goFlux]{empirical.prec}}
#'
#' @examples
#' data(manID.UGGA)
#' CH4_flux <- goFlux(manID.UGGA, "CH4dry_ppb")
#' CH4_best <- best.flux(CH4_flux)
#' CH4_qc <- qc.flags(CH4_best, dataframe = manID.UGGA, gastype = "CH4dry_ppb")
#' CH4_qc[, c("UniqueID", "best.flux", "qc.c0", "qc.convex", "qc.min.obs",
#'            "qc.ambient", "qc.ambient.dev", "qc.clock", "qc.noisy", "qc.any")]
#' @export
qc.flags <- function(flux.result, dataframe = NULL, gastype = NULL,
                     c0.mult = 1.5, by = NULL, min.obs = 60,
                     ambient.sigma = 3, ambient.secs = 10, ambient.pre = 60,
                     seal.time = "start.time", noisy.mult = 1.5) {

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
  for(a in c("ambient.secs", "ambient.pre")){
    v <- get(a)
    if(!is.numeric(v) || length(v) != 1 || !(v > 0)) stop("'", a, "' must be a single number > 0")}
  if(!is.character(seal.time) || length(seal.time) != 1) stop("'seal.time' must be a character string")

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
  if(!is.null(dataframe) && (!is.null(ambient.sigma) || !is.null(noisy.mult))){
    dd <- as.data.frame(dataframe)
    dd$qc_by <- g[match(as.character(dd$UniqueID), uid)]
    dd <- dd[!is.na(dd$qc_by), ]
    # per-closure and group second-difference precision (empirical.prec);
    # each closure gets its group's value at the closure's own logging interval
    ep <- empirical.prec(dd, gastype, by = "qc_by")
    clo <- attr(ep, "closures")
    i.c <- match(uid, clo$UniqueID)
    s.clo <- clo$prec[i.c]; dt.clo <- clo$dt_s[i.c]
    s.grp <- vapply(seq_len(n), function(i){
      e <- ep[ep$qc_by == g[i], , drop = FALSE]
      if(nrow(e) == 0) return(NA_real_)
      j <- if(is.na(dt.clo[i])) integer(0) else which(!is.na(e$dt_s) & e$dt_s == dt.clo[i])
      if(length(j) == 0) j <- which.max(ifelse(is.na(e$n.closures), -1, e$n.closures))
      e$prec[j[1]] }, numeric(1))
    fx$qc.prec <- s.grp

    # qc.ambient (at sealing) and qc.clock
    has.seal <- any(grepl(paste("\\<", seal.time, "\\>", sep = ""), names(dd))) &&
      any(grepl("\\<POSIX.time\\>", names(dd)))
    if(!is.null(ambient.sigma) && has.seal){
      sp <- split(dd, as.character(dd$UniqueID))
      res <- t(vapply(seq_len(n), function(i){
        x <- sp[[uid[i]]]
        if(is.null(x)) return(c(NA_real_, NA_real_, NA_real_))
        tm <- as.numeric(x$POSIX.time); v <- x[[gastype]]
        seal <- as.numeric(x[[seal.time]]); seal <- seal[!is.na(seal)][1]
        if(is.na(seal)) return(c(NA_real_, NA_real_, NA_real_))
        w <- tm[!is.na(x$flag) & x$flag == 1]
        dt <- if(is.finite(dt.clo[i])) dt.clo[i] else 1
        clock <- if(length(w)) as.numeric(min(w) < seal - dt) else NA_real_
        if(isTRUE(clock == 1)) return(c(NA_real_, NA_real_, 1))
        pre <- is.finite(v) & tm >= seal - ambient.pre & tm < seal
        sl <- is.finite(v) & tm >= seal & tm < seal + ambient.secs
        if(sum(pre) < 3 || sum(sl) < 1) return(c(NA_real_, NA_real_, clock))
        amb <- stats::median(v[pre])
        tol <- max(s.grp[i], stats::mad(v[pre], constant = 1.4826), na.rm = TRUE)
        if(!is.finite(tol) || tol <= 0) return(c(NA_real_, NA_real_, clock))
        tp <- tm[pre] - seal
        drift <- unname(stats::coef(stats::lm(v[pre] ~ tp))[2]) * ambient.pre
        c((stats::median(v[sl]) - amb) / tol, drift / tol, clock) }, numeric(3)))
      fx$qc.ambient.dev <- res[, 1]
      fx$qc.ambient.drift <- res[, 2]
      fx$qc.ambient <- ifelse(is.na(res[, 1]), NA, abs(res[, 1]) > ambient.sigma)
      fx$qc.clock <- ifelse(is.na(res[, 3]), NA, res[, 3] == 1)
      flags <- c(flags, "qc.ambient", "qc.clock")
    }

    # qc.noisy: closure precision vs group precision
    if(!is.null(noisy.mult)){
      fx$qc.noisy.ratio <- s.clo / s.grp
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
