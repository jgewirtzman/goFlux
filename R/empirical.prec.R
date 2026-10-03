#' Empirical instrument precision from the closures themselves
#'
#' Estimates the instrument precision (1 sigma) of a gas measurement from the
#' data, instead of the manufacturer specification. The default estimator
#' works closure by closure on the second differences of the gas inside the
#' fitted window, which cancel the flux (any straight line) exactly, and
#' summarises a group of closures (an analyzer on a field day, say) by the
#' median of the per-closure values.
#'
#' @param dataframe data.frame; output from the ID functions
#'                  (\code{\link[goFlux]{obs.win}} followed by
#'                  \code{\link[goFlux]{click.peak2}} or
#'                  \code{\link[goFlux]{crop.meas}}). Must contain
#'                  \code{gastype}; for \code{method = "hadamard"} and
#'                  \code{"allan"} it must also contain \code{UniqueID} and
#'                  \code{flag}. \code{method = "mad_diff1"} also accepts a
#'                  whole imported record without \code{UniqueID}.
#' @param gastype character string; the gas column, e.g. \code{"CH4dry_ppb"}.
#' @param method character string; \code{"hadamard"} (default), the
#'               per-closure second-difference precision summarised by its
#'               median over the closures of each group (see Details);
#'               \code{"allan"}, \eqn{SD(\Delta x) / \sqrt{2}} per
#'               \code{UniqueID} on the rows with \code{flag == 1} (the Allan
#'               deviation at one sampling interval); \code{"mad_diff1"},
#'               superseded, \eqn{MAD(\Delta x) \times 1.4826 / \sqrt{2}}
#'               over the whole record (or \code{by} group), kept for
#'               comparison with earlier results only (see Details).
#'               \code{method = "mad"}, the name used before, partially
#'               matches \code{"mad_diff1"}.
#' @param by character string; an optional column of \code{dataframe} defining
#'           the groups (e.g. analyzer x field day) over which the precision
#'           is summarised. Default \code{NULL}: one group, all closures.
#'           Ignored for \code{method = "allan"}, which is always per
#'           \code{UniqueID}. Use \code{by = "UniqueID"} for one row per
#'           closure (the per-closure values are also returned as the
#'           attribute \code{"closures"}).
#' @param tol numerical; the fractional change in the spacing of
#'            \code{POSIX.time} that starts a new run of constant logging
#'            interval (default 0.2, i.e. 20 \%). Differences never span two
#'            runs. See Details.
#' @param warn logical; give the diagnostic warnings described under
#'             \emph{Checks}? Default \code{TRUE}. The check values are
#'             returned as columns either way.
#'
#' @details
#' \strong{Estimator (\code{method = "hadamard"}).} For each closure
#' (\code{UniqueID}) the rows with \code{flag == 1} are put in time order and
#' the precision is
#' \deqn{\sigma_i = MAD(x_{k+1} - 2 x_k + x_{k-1}) / \sqrt{6},}{sigma_i = MAD(x[k+1] - 2 x[k] + x[k-1]) / sqrt(6),}
#' where MAD is the normal-consistent median absolute deviation about the
#' median (\code{stats::mad}, constant 1.4826). A straight line, i.e. the flux,
#' cancels exactly in a second difference, and slow curvature leaves only a
#' term of order curvature x \eqn{\Delta t^2}; white noise of standard
#' deviation \eqn{\sigma} gives \eqn{Var = 6 \sigma^2}. This is the
#' one-sample Hadamard deviation, the drift-insensitive counterpart of the
#' Allan deviation (Riley 2008). A closure needs at least 6 second differences;
#' otherwise its \eqn{\sigma_i} is \code{NA}. The precision of a group
#' (\code{by}) is the \emph{median} of \eqn{\sigma_i} over its closures, so that
#' a minority of very large, curved or bubbly closures cannot set the group
#' value.
#'
#' \strong{Why not the whole-record first-difference MAD.} The superseded
#' \code{method = "mad_diff1"} pools the first differences of a whole record
#' and centres them on a single median. A first difference contains the
#' closure's slope x \eqn{\Delta t}, and closures differ in slope, so the
#' spread of slopes between closures is counted as instrument noise. Two
#' noise-free closures rising 0 and 10 ppb per step give a "precision" of
#' about 7 ppb. The inflation is small at 1 Hz but reaches 1.4-2.3x at 5-10 s
#' logging, where the slope part of a difference is 5-10 times larger
#' (\url{https://github.com/jgewirtzman/fluxqc/issues/1}). It is kept for
#' comparison with earlier results only.
#'
#' \strong{Logging interval.} A second difference must not span a change of
#' logging interval or a gap. When \code{dataframe} contains
#' \code{POSIX.time}, each closure is split into runs of constant logging
#' interval: a run breaks where the spacing changes by more than \code{tol};
#' duplicate or backwards timestamps break a run; differences are taken within
#' runs only. A closure's interval \code{dt_s} is the median spacing of its
#' flagged rows (nearest 0.5 s), and only differences at that interval are
#' used. If the closures of a group were logged at different intervals, one
#' row per interval is returned and a warning is given: pick the row whose
#' \code{dt_s} matches the closures you are evaluating. Without
#' \code{POSIX.time} the rows of each closure are taken as one run in their
#' present order and \code{dt_s} is \code{NA}. For \code{"mad_diff1"} the same
#' run logic applies to the whole record (runs shorter than 3 differences,
#' and intervals seen fewer than 3 times, are dropped).
#'
#' \strong{Checks} (returned as columns; \code{warn = TRUE} also warns):
#' \describe{
#'   \item{\code{ac1}}{lag-1 autocorrelation of the second differences
#'     (median over closures). For white noise it is -2/3. Above -0.5 the noise
#'     is red (random-walk or drift-dominated): neighbouring samples are not
#'     independent, the second differences understate the noise that matters
#'     over a closure, and an Allan plot is needed.}
#'   \item{\code{zero.frac}}{fraction of first differences that are exactly
#'     zero (pooled over closures). At 0.3 or more the analyzer probably logs
#'     one gas per row and carries the other forward (e.g. the Picarro G4301,
#'     whose CH4 is fresh on every other row only), or the record is quantized
#'     at a resolution close to the noise. Either way a difference-based
#'     precision is unreliable: first differences are biased low (by 10-40x
#'     for alternate-row loggers), and second differences no longer cancel the
#'     slope, because the fresh values are two rows apart. Keep each gas's
#'     fresh rows before estimating, or quote the resolution.}
#'   \item{\code{prec.d1c}, \code{d1c.ratio}}{the precision from first
#'     differences centred on each closure's own median
#'     (\eqn{MAD / \sqrt{2}}{MAD / sqrt(2)} of the pooled centred differences)
#'     and its ratio to \code{prec}. For white noise around straight lines the
#'     ratio is 1; above 1.2 (warned) curvature or trend leaks into the first
#'     differences, or the noise is not white; below 0.8 (warned) points at
#'     alternate-row or quantized data.}
#' }
#'
#' \strong{Drift is a separate problem.} Any difference-based precision,
#' this one included, is blind to instrument drift and to slow leaks: both
#' change the concentration smoothly over a closure and so look like a flux,
#' not like noise. Bound them with periodic blank closures on an inert surface
#' (a sealed chamber on a plate or foil), whose apparent flux is the drift or
#' leak floor; the trend of the ambient record before and after closures can
#' flag an unstable analyzer but cannot measure drift inside the chamber.
#'
#' \strong{Datasheet comparison.} Manufacturers quote precision at 1 s. A record
#' logged every \eqn{\Delta t} seconds averages white noise down by
#' \eqn{\sqrt{\Delta t}}{sqrt(dt)}, so before comparing the \code{prec} of a
#' 5 s or 10 s record with the datasheet, divide the datasheet value by
#' \eqn{\sqrt{\Delta t}}{sqrt(dt)} (\code{\link[goFlux]{spec.at.interval}};
#' \code{dt_s} is returned for this). Otherwise an analyzer can appear to beat
#' its datasheet when it is several times worse. The rescaling assumes white
#' noise and fails beyond the analyzer's Allan minimum.
#'
#' \strong{Use in the MDF.} Pass the result to the \code{prec} argument of the
#' import functions (e.g. \code{\link[goFlux]{import.LI7810}}) or write it into
#' the \code{*_prec} columns of the data frame before
#' \code{\link[goFlux]{goFlux}}: the \code{prec} column of the output and
#' \code{\link[goFlux]{MDF}} then use it, MDF = z x prec / t x flux.term with
#' \code{t} in seconds. See \code{\link[goFlux]{MDF}} for the confidence level
#' (\code{conf}).
#'
#' @returns A data.frame. For \code{method = "hadamard"}: one row per group (and
#'          per logging interval, see Details) with the grouping column (if
#'          any), \code{prec} (median of the per-closure precisions, units of
#'          \code{gastype}), \code{n.closures} (closures with a finite
#'          \eqn{\sigma_i}), \code{dt_s} (logging interval, s),
#'          \code{n.diff} (second differences used), \code{ac1},
#'          \code{zero.frac}, \code{prec.d1c} and \code{d1c.ratio} (see
#'          \emph{Checks}). The per-closure values (\code{UniqueID}, the
#'          grouping column, \code{dt_s}, \code{prec}, \code{n.diff},
#'          \code{ac1}, \code{zero.frac}, \code{prec.d1c}, \code{d1c.ratio})
#'          are in the attribute \code{"closures"}. For \code{"allan"}: one row
#'          per \code{UniqueID} with \code{prec} and \code{n} (differences
#'          used). For \code{"mad_diff1"}: one row per record or group and
#'          logging interval with \code{prec}, \code{n} (differences used) and
#'          \code{dt_s}.
#'
#' @references
#' Riley, W. J. (2008). Handbook of Frequency Stability Analysis. NIST Special
#' Publication 1065. National Institute of Standards and Technology.
#' \doi{10.6028/NIST.SP.1065}
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{MDF}}, \code{\link[goFlux]{goFlux}},
#'          \code{\link[goFlux]{best.flux}}, \code{\link[goFlux]{spec.at.interval}}
#'
#' @examples
#' data(manID.UGGA)
#' # second-difference precision of the closure (one group), with the checks
#' empirical.prec(manID.UGGA, "CH4dry_ppb")
#' # per closure
#' attr(empirical.prec(manID.UGGA, "CH4dry_ppb"), "closures")
#' # Allan deviation on the flagged window
#' empirical.prec(manID.UGGA, "CH4dry_ppb", method = "allan")
#' # two closures logged every 10 s with very different slopes: the median
#' # second-difference precision recovers the noise (1 ppb); the superseded
#' # whole-record first-difference MAD counts the slope difference as noise
#' set.seed(1)
#' t0 <- as.POSIXct("2024-06-01", tz = "UTC")
#' two <- data.frame(UniqueID = rep(c("slow", "fast"), each = 60), flag = 1,
#'                   POSIX.time = t0 + 10 * (0:119),
#'                   CH4dry_ppb = 2000 + c(0 * (1:60), 10 * (1:60)) + rnorm(120, 0, 1))
#' empirical.prec(two, "CH4dry_ppb")$prec
#' empirical.prec(two, "CH4dry_ppb", method = "mad_diff1")$prec
#' @export
empirical.prec <- function(dataframe, gastype, method = c("hadamard", "allan", "mad_diff1"),
                           by = NULL, tol = 0.2, warn = TRUE) {

  # Check arguments
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.data.frame(dataframe)) stop("'dataframe' must be of class data.frame")
  if(missing(gastype)) stop("'gastype' is required")
  if(!is.character(gastype) | length(gastype) != 1) stop("'gastype' must be a character string")
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe)))){
    stop("'dataframe' must contain a column that matches 'gastype'")}
  method <- match.arg(method)
  if(!is.null(by)){
    if(!is.character(by) | length(by) != 1) stop("'by' must be a character string")
    if(!any(grepl(paste("\\<", by, "\\>", sep = ""), names(dataframe)))){
      stop("'dataframe' must contain a column that matches 'by'")}}
  if(!is.numeric(tol) | length(tol) != 1) stop("'tol' must be a single number")
  if(tol < 0) stop("'tol' must be >= 0")
  if(!is.logical(warn) | length(warn) != 1) stop("'warn' must be TRUE or FALSE")

  x <- dataframe[[gastype]]
  has.time <- any(grepl("\\<POSIX.time\\>", names(dataframe)))
  p <- if(has.time) as.numeric(dataframe$POSIX.time) else NULL

  if(method == "hadamard"){
    for(col in c("UniqueID", "flag")){
      if(!any(grepl(paste("\\<", col, "\\>", sep = ""), names(dataframe)))){
        stop("'dataframe' must contain the column '", col, "' for method = 'hadamard' ",
             "(the precision is estimated inside the fitted windows)")}}
    keep <- !is.na(dataframe$flag) & dataframe$flag == 1
    if(!any(keep)) stop("'dataframe' has no rows with flag == 1")
    uid <- as.character(dataframe$UniqueID)[keep]
    xv <- x[keep]; tv <- if(has.time) p[keep] else NULL
    gv <- if(is.null(by)) rep("all", length(uid)) else as.character(dataframe[[by]])[keep]
    idx <- split(seq_along(uid), uid)
    cl <- lapply(idx, function(i) hadamard.closure(xv[i], if(has.time) tv[i] else NULL, tol))
    num <- function(nm) vapply(cl, function(r) as.numeric(r[[nm]]), numeric(1), USE.NAMES = FALSE)
    clo <- data.frame(UniqueID = names(idx),
                      group = vapply(idx, function(i) gv[i[1]], "", USE.NAMES = FALSE),
                      dt_s = num("dt_s"), prec = num("prec"), n.diff = as.integer(num("n2")),
                      ac1 = num("ac1"), zero.frac = num("zero.frac"), prec.d1c = num("prec.d1c"),
                      stringsAsFactors = FALSE)
    clo$d1c.ratio <- clo$prec.d1c / clo$prec
    # group summary, per group and logging interval
    key <- paste(clo$group, ifelse(is.na(clo$dt_s), "NA", clo$dt_s), sep = "\r")
    ks <- unique(key)
    ks <- ks[order(clo$group[match(ks, key)], clo$dt_s[match(ks, key)])]
    out <- do.call(rbind, lapply(ks, function(k){
      j <- which(key == k)
      d1 <- unlist(lapply(cl[j], function(r) r$d1c))
      nz <- sum(num("n.zero")[j]); nd <- sum(num("n1")[j])
      s <- if(any(is.finite(clo$prec[j]))) stats::median(clo$prec[j], na.rm = TRUE) else NA_real_
      a <- if(any(is.finite(clo$ac1[j]))) stats::median(clo$ac1[j], na.rm = TRUE) else NA_real_
      s.d1 <- if(length(d1) >= 3) stats::mad(d1, constant = 1.4826) / sqrt(2) else NA_real_
      data.frame(group = clo$group[j[1]], prec = s, n.closures = sum(is.finite(clo$prec[j])),
                 dt_s = clo$dt_s[j[1]], n.diff = sum(clo$n.diff[j]), ac1 = a,
                 zero.frac = if(nd > 0) nz / nd else NA_real_, prec.d1c = s.d1,
                 d1c.ratio = s.d1 / s, stringsAsFactors = FALSE) }))
    rownames(out) <- NULL
    if(warn){
      lab <- function(i) if(is.null(by)) "" else paste0("group '", out$group[i], "': ")
      for(g in unique(out$group)){
        r <- out[out$group == g, ]
        if(nrow(r) > 1){
          warning(if(is.null(by)) "closures have " else paste0("group '", g, "' has "),
                  nrow(r), " logging intervals (", paste(r$dt_s, collapse = ", "),
                  " s); precision is reported per interval (one row each)", call. = FALSE)}}
      for(i in seq_len(nrow(out))){
        if(!is.na(out$zero.frac[i]) && out$zero.frac[i] >= 0.3){
          warning(lab(i), round(100 * out$zero.frac[i]), " % of first differences of ",
                  gastype, " are exactly zero: the analyzer probably logs one gas per row ",
                  "(alternate-row logging, e.g. Picarro G4301) or the record is quantized; ",
                  "difference-based precision is unreliable. Keep each gas's fresh rows first.", call. = FALSE)}
        if(!is.na(out$ac1[i]) && out$ac1[i] > -0.5){
          warning(lab(i), "lag-1 autocorrelation of second differences is ", round(out$ac1[i], 2),
                  " (white noise: -0.67): red or drift-dominated noise; the second-difference ",
                  "precision understates the noise over a closure (check an Allan plot).", call. = FALSE)}
        if(!is.na(out$d1c.ratio[i]) && (out$d1c.ratio[i] > 1.2 || out$d1c.ratio[i] < 0.8)){
          warning(lab(i), "precision from centred first differences is ", round(out$d1c.ratio[i], 2),
                  "x the second-difference precision (expected ~1): ",
                  if(out$d1c.ratio[i] > 1.2) "trend or curvature leaks into first differences, or the noise is not white."
                  else "alternate-row or quantized data?", call. = FALSE)}
      }
    }
    if(is.null(by)){ out$group <- NULL; clo$group <- NULL }
    else {
      names(out)[names(out) == "group"] <- by
      if(by == "UniqueID") clo$group <- NULL else names(clo)[names(clo) == "group"] <- by }
    attr(out, "closures") <- clo
  } else if(method == "mad_diff1"){
    # Superseded: MAD of first differences per run of constant logging
    # interval over the whole record. Returns one row per interval: prec, n
    # (differences used) and dt_s. Without timestamps: one pooled MAD.
    f.mad <- function(v, tm){
      if(is.null(tm)){
        v <- v[!is.na(v)]; d <- diff(v)
        if(length(d) < 2) return(data.frame(prec = NA_real_, n = length(d), dt_s = NA_real_))
        return(data.frame(prec = stats::mad(d, constant = 1.4826) / sqrt(2),
                          n = length(d), dt_s = NA_real_))
      }
      ok <- is.finite(v) & !is.na(tm); v <- v[ok]; tm <- tm[ok]
      d <- diff(tm); dx <- diff(v)
      # duplicate or backwards timestamps are not logging intervals
      valid <- is.finite(d) & d > 0
      none <- data.frame(prec = NA_real_, n = 0L, dt_s = NA_real_)
      if(sum(valid) < 3) return(none)
      dv <- d[valid]
      # a run breaks where the spacing changes by more than tol (fraction)
      brk <- c(TRUE, abs(diff(dv)) > tol * pmax(dv[-1], 1e-9))
      run <- cumsum(brk)
      # a lone gap between files is not a run
      keep <- run %in% which(tabulate(run) >= 3)
      idx <- which(valid)[keep]
      if(length(idx) < 3) return(none)
      dt_r <- round(d[idx] * 2) / 2   # nearest 0.5 s
      # one pooled MAD per interval (runs with the same interval are pooled)
      out <- do.call(rbind, lapply(split(idx, dt_r), function(i)
        data.frame(prec = stats::mad(dx[i], constant = 1.4826) / sqrt(2),
                   n = length(i), dt_s = round(d[i[1]] * 2) / 2)))
      # an interval with fewer than 3 differences is not a logging interval
      out <- out[out$n >= 3, , drop = FALSE]
      if(nrow(out) == 0) return(none)
      rownames(out) <- NULL
      out
    }
    if(is.null(by)){
      out <- f.mad(x, p)
      if(nrow(out) > 1 && warn){
        warning("record has ", nrow(out), " logging intervals (",
                paste(out$dt_s, collapse = ", "), " s); precision is reported ",
                "per interval (one row each)", call. = FALSE)}
    } else {
      g <- as.character(dataframe[[by]])
      idx.g <- split(seq_along(x), g)
      r <- lapply(names(idx.g), function(k){
        i <- idx.g[[k]]
        o <- f.mad(x[i], if(is.null(p)) NULL else p[i])
        if(nrow(o) > 1 && warn){
          warning("group '", k, "' has ", nrow(o), " logging intervals (",
                  paste(o$dt_s, collapse = ", "), " s); precision is reported ",
                  "per interval (one row each)", call. = FALSE)}
        cbind(g = rep(k, nrow(o)), o, stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, r)
      rownames(out) <- NULL
      names(out)[1] <- by
    }
  } else {
    if(!any(grepl("\\<UniqueID\\>", names(dataframe)))){
      stop("'dataframe' must contain the column 'UniqueID' for method = 'allan'")}
    if(!any(grepl("\\<flag\\>", names(dataframe)))){
      stop("'dataframe' must contain the column 'flag' for method = 'allan'")}
    d <- dataframe[!is.na(dataframe$flag) & dataframe$flag == 1, ]
    f.allan <- function(v){ v <- v[!is.na(v)]; dd <- diff(v)
      if(length(dd) < 2) return(c(NA_real_, length(dd)))
      c(stats::sd(dd) / sqrt(2), length(dd)) }
    r <- t(sapply(split(d[[gastype]], as.character(d$UniqueID)), f.allan))
    out <- data.frame(UniqueID = rownames(r), prec = r[, 1], n = as.integer(r[, 2]),
                      row.names = NULL, stringsAsFactors = FALSE)
  }
  out
}

# Second-difference (one-sample Hadamard) precision of ONE closure, with the
# checks. v: gas values of the flagged rows; tm: their POSIX.time as numeric
# (NULL: rows taken in their order as one run); tol: fractional spacing change
# that breaks a run. Differences never span a change of interval, a gap or a
# duplicate/backwards timestamp, and only differences at the closure's own
# interval (median spacing) are used.
hadamard.closure <- function(v, tm = NULL, tol = 0.2){
  na <- list(prec = NA_real_, n2 = 0L, dt_s = NA_real_, ac1 = NA_real_, zero.frac = NA_real_,
             prec.d1c = NA_real_, d1c = numeric(0), n.zero = 0, n1 = 0)
  if(is.null(tm)){
    v <- v[is.finite(v)]
    if(length(v) < 3) return(na)
    dx <- diff(v); run <- rep(1L, length(dx)); dt_s <- NA_real_
  } else {
    ok <- is.finite(v) & is.finite(tm); v <- v[ok]; tm <- tm[ok]
    o <- order(tm); v <- v[o]; tm <- tm[o]
    if(length(v) < 3) return(na)
    d <- diff(tm); dx <- diff(v)
    valid <- d > 0
    if(!any(valid)) return(na)
    dm <- stats::median(d[valid])
    dt_s <- if(dm >= 0.5) round(dm * 2) / 2 else signif(dm, 2)
    # a run breaks after an invalid spacing or where the spacing changes by > tol
    brk <- c(TRUE, abs(diff(d)) > tol * pmax(d[-1], 1e-9) | !valid[-length(valid)])
    run <- cumsum(brk)
    # only differences at the closure's own interval
    run[!(valid & abs(d - dm) <= tol * dm)] <- NA_integer_
  }
  d1 <- dx[!is.na(run)]
  n1 <- length(d1)
  k <- length(dx)
  same <- if(k >= 2) !is.na(run[-1]) & !is.na(run[-k]) & run[-1] == run[-k] else logical(0)
  d2 <- if(k >= 2) (dx[-1] - dx[-k])[same] else numeric(0)
  pos <- which(same)                      # pair (pos, pos + 1) of first differences
  r2 <- if(k >= 2) run[-1][same] else integer(0)
  n2 <- length(d2)
  prec <- if(n2 >= 6) stats::mad(d2, constant = 1.4826) / sqrt(6) else NA_real_
  # lag-1 autocorrelation of consecutive second differences within a run
  ac1 <- NA_real_
  if(n2 >= 8){
    adj <- which(diff(pos) == 1 & r2[-1] == r2[-n2])
    if(length(adj) >= 5){
      a <- d2[adj]; b <- d2[adj + 1]
      if(stats::sd(a) > 0 && stats::sd(b) > 0) ac1 <- stats::cor(a, b)
    }
  }
  # first differences centred on the closure's own median
  d1c <- if(n1 >= 3) d1 - stats::median(d1) else numeric(0)
  prec.d1c <- if(n1 >= 3) stats::mad(d1, constant = 1.4826) / sqrt(2) else NA_real_
  list(prec = prec, n2 = as.integer(n2), dt_s = dt_s, ac1 = ac1,
       zero.frac = if(n1 > 0) sum(d1 == 0) / n1 else NA_real_,
       prec.d1c = prec.d1c, d1c = d1c, n.zero = sum(d1 == 0), n1 = n1)
}
