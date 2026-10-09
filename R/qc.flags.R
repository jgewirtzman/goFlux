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
#' @param min.secs numerical; \code{qc.min.secs} fires when the closure
#'                duration (\code{\link[goFlux]{closure.time}} of the rows with
#'                \code{flag == 1} in \code{dataframe}) is below it, in
#'                seconds, whatever the logging interval. Default \code{NULL}
#'                (skipped); needs \code{dataframe}.
#' @param convex.p numerical; significance level of the curvature test of
#'                \code{qc.convex} (see Details). Default 0.05.
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
#'                \code{NULL}: \code{"cham.close"} when \code{dataframe} has
#'                it (output of \code{\link[goFlux]{crop.meas}} or
#'                \code{\link[goFlux]{windows.from.table}}, whose
#'                \code{start.time} has been moved to the start of the window),
#'                else \code{"start.time"} (the auxfile start carried by
#'                \code{\link[goFlux]{obs.win}} and
#'                \code{\link[goFlux]{click.peak2}}). Never use the start of
#'                the fitting window (\code{start.time_corr}).
#' @param noisy.mult numerical; \code{qc.noisy} fires when the closure's own
#'                precision exceeds \code{noisy.mult} times the group precision
#'                (both from \code{\link[goFlux]{empirical.prec}}, second
#'                differences of the flagged rows; the group value is the median
#'                over its closures, at the closure's own logging interval): a
#'                disturbed closure. Default 1.5. Needs \code{dataframe}.
#' @param leak.rate numerical; optional bench-measured leak rate, the fraction
#'                of the headspace-ambient difference exchanged per second
#'                (e.g. the decay constant of a spiked, sealed chamber).
#'                Default \code{NULL}. See \emph{Leaks}.
#' @param blank.slope numerical; optional apparent slope (units of
#'                \code{gastype} per second) of blank closures on an inert
#'                surface, which measures drift and leaks together. Default
#'                \code{NULL}. See \emph{Leaks}.
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
#' fluxes with (\code{blank.slope}); check the analyzer with reference, zero
#' or span gases. The ambient shoulders before a closure mostly reflect real
#' ambient variability: their trend (\code{qc.ambient.drift}) is a flag only,
#' it cannot measure drift inside the chamber and must not be subtracted.
#'
#' \strong{Leaks (\code{qc.leak}).} With \code{leak.rate}, the apparent flux a
#' leak alone could produce is \code{leak.rate} times the largest
#' headspace-ambient difference during the closure (the concentration range of
#' the flagged rows of \code{dataframe}, else \code{|Ct - C0|}) times
#' \code{flux.term}; with \code{blank.slope}, \code{|blank.slope|} times
#' \code{flux.term}; the larger of the two when both are given
#' (\code{qc.leak.flux}). \code{qc.leak} fires when it reaches the minimal
#' detectable flux of the row (\code{det.MDF} from
#' \code{\link[goFlux]{flux.class}} if present, else goFlux's \code{MDF}): the
#' flux cannot be separated from leakage or drift. Skipped when both are
#' \code{NULL}.
#'
#' \strong{Convex trace (\code{qc.convex}).} An accelerating concentration
#' change, which no chamber mechanism produces, points at a leak, a disturbance
#' or a delayed seal. With \code{dataframe}, a quadratic is fitted to the rows
#' with \code{flag == 1} of each closure (concentration on \code{Etime}); the
#' flag fires when the quadratic term has the sign of the overall linear slope
#' and its p-value is below \code{convex.p} (at least 6 rows; \code{NA}
#' otherwise). Without \code{dataframe} it falls back to the
#' Hutchinson-Mosier curvature, \code{HM.k < 0}, which can only occur when
#' \code{\link[goFlux]{goFlux}} was run with \code{k.min < 0}.
#'
#' The CO2 tracer test (CO2 must accumulate in a sealed chamber on a respiring
#' surface) is not a flag here because it is a plain comparison between two
#' \code{best.flux} outputs: see \code{\link[goFlux]{co2.tracer}}.
#'
#' @returns \code{flux.result} with the columns \code{qc.c0}, \code{qc.convex},
#'          \code{qc.min.obs}, \code{qc.min.secs}, \code{qc.ambient}, \code{qc.clock},
#'          \code{qc.noisy}, \code{qc.leak} (those requested) and \code{qc.any}
#'          appended, plus
#'          the helper values \code{qc.c0.ratio}, \code{qc.ambient.dev} and
#'          \code{qc.ambient.drift} (signed, in units of the ambient
#'          tolerance), \code{qc.prec} (the group precision used) and
#'          \code{qc.noisy.ratio} and \code{qc.leak.flux}.
#'
#' @include goFlux-package.R
#' @include empirical.prec.R
#' @include closure.time.R
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
                     c0.mult = 1.5, by = NULL, min.obs = 60, min.secs = NULL,
                     convex.p = 0.05,
                     ambient.sigma = 3, ambient.secs = 10, ambient.pre = 60,
                     seal.time = NULL, noisy.mult = 1.5,
                     leak.rate = NULL, blank.slope = NULL) {

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
  if(is.null(seal.time)){
    seal.time <- if(!is.null(dataframe) && "cham.close" %in% names(dataframe)) "cham.close" else "start.time"}
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

  # flagged rows of each closure, in time order (for qc.convex and qc.min.secs)
  win <- NULL
  if(!is.null(dataframe)){
    dw <- as.data.frame(dataframe)
    dw <- dw[!is.na(dw$flag) & dw$flag == 1, ]
    dw <- dw[order(dw$UniqueID, dw$Etime), ]
    win <- split(dw, as.character(dw$UniqueID))
  }

  # qc.convex: accelerating change (quadratic term with the sign of the slope)
  if(!is.null(win)){
    cv <- vapply(uid, function(u){
      d <- win[[u]]
      if(is.null(d) || nrow(d) < 6) return(NA)
      y <- d[[gastype]]; tt <- d$Etime
      m <- try(stats::lm(y ~ tt + I(tt^2)), silent = TRUE)
      if(inherits(m, "try-error")) return(NA)
      co <- suppressWarnings(summary(m))$coefficients   # a noise-free trace warns
      if(nrow(co) < 3) return(NA)
      net <- unname(stats::coef(stats::lm(y ~ tt))[2])
      sign(co[3, 1]) == sign(net) && co[3, 4] < convex.p }, logical(1), USE.NAMES = FALSE)
    fx$qc.convex <- cv
  } else {
    fx$qc.convex <- if(any(grepl("\\<HM.k\\>", names(fx)))) ifelse(is.na(fx$HM.k), NA, fx$HM.k < 0) else NA
  }
  flags <- c(flags, "qc.convex")

  # qc.min.secs: closure duration in seconds
  if(!is.null(min.secs)){
    if(is.null(win)) stop("'min.secs' needs 'dataframe'")
    dur <- vapply(uid, function(u){
      d <- win[[u]]
      if(is.null(d) || nrow(d) == 0) return(NA_real_)
      closure.time(d$Etime) }, numeric(1), USE.NAMES = FALSE)
    fx$qc.min.secs <- ifelse(is.na(dur), NA, dur < min.secs)
    flags <- c(flags, "qc.min.secs")
  }

  # qc.min.obs
  if(!is.null(min.obs)){
    fx$qc.min.obs <- if(any(grepl("\\<nb.obs\\>", names(fx)))) ifelse(is.na(fx$nb.obs), NA, fx$nb.obs < min.obs) else NA
    flags <- c(flags, "qc.min.obs")
  }

  # flags needing the concentration data
  if(!is.null(dataframe) && (!is.null(ambient.sigma) || !is.null(noisy.mult))){
    dd <- as.data.frame(dataframe)
    # per-closure and group second-difference precision (empirical.prec);
    # each closure gets its group's value at the closure's own logging interval
    pc <- prec.by.closure(dd, gastype, uid, g)
    s.clo <- pc$closure; s.grp <- pc$group; dt.clo <- pc$dt
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

  # qc.leak: apparent flux from a measured leak rate or blank-closure slope
  if(!is.null(leak.rate) || !is.null(blank.slope)){
    if(!any(grepl("\\<flux.term\\>", names(fx)))) stop("'flux.result' must contain 'flux.term' for the leak flag")
    app <- rep(0, n)
    if(!is.null(leak.rate)){
      if(!is.numeric(leak.rate) || length(leak.rate) != 1 || leak.rate < 0) stop("'leak.rate' must be a single number >= 0")
      dC <- if(!is.null(dataframe)){
        dd <- as.data.frame(dataframe)
        dd <- dd[!is.na(dd$flag) & dd$flag == 1, ]
        r <- vapply(split(dd[[gastype]], as.character(dd$UniqueID)),
                    function(v) diff(range(v, na.rm = TRUE)), numeric(1))
        r[match(uid, names(r))]
      } else if(all(c("C0", "Ct") %in% names(fx))) abs(fx$Ct - fx$C0) else rep(NA_real_, n)
      app <- leak.rate * dC * fx$flux.term
    }
    if(!is.null(blank.slope)){
      if(!is.numeric(blank.slope) || length(blank.slope) != 1) stop("'blank.slope' must be a single number")
      app <- pmax(app, abs(blank.slope) * fx$flux.term)
    }
    mdf <- if("det.MDF" %in% names(fx)) fx$det.MDF else if("MDF" %in% names(fx)) fx$MDF else rep(NA_real_, n)
    fx$qc.leak.flux <- app
    fx$qc.leak <- ifelse(is.na(app) | is.na(mdf), NA, app >= mdf)
    flags <- c(flags, "qc.leak")
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
