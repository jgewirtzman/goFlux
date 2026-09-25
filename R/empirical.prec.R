#' Empirical instrument precision from first differences
#'
#' Estimates the instrument precision (1 sigma) of a gas measurement from the
#' data themselves, instead of the manufacturer specification. First
#' differences remove any smooth trend, so the flux signal cancels and only the
#' instrument noise remains: for white noise \eqn{Var(x_{i+1} - x_i) = 2\sigma^2}.
#'
#' @param dataframe data.frame; output from import or ID functions
#'                  (\code{\link[goFlux]{obs.win}}, \code{\link[goFlux]{click.peak2}},
#'                  \code{\link[goFlux]{crop.meas}}) or a whole imported record.
#'                  Must contain \code{gastype}; for \code{method = "allan"} it
#'                  must also contain \code{UniqueID} and \code{flag}.
#' @param gastype character string; the gas column, e.g. \code{"CH4dry_ppb"}.
#' @param method character string; \code{"mad"} (default) is
#'               \eqn{MAD(\Delta x) \times 1.4826 / \sqrt{2}} over the whole
#'               record (or per group, see \code{by}): robust to the isolated
#'               large differences produced by moves, breaths and chamber
#'               changes, so the entire field record can be used without
#'               selecting quiet periods. \code{"allan"} is
#'               \eqn{SD(\Delta x) / \sqrt{2}}, the Allan deviation at one
#'               sampling interval, computed per \code{UniqueID} on the rows
#'               with \code{flag == 1}.
#' @param by character string; for \code{method = "mad"}, an optional column of
#'           \code{dataframe} defining the field periods (e.g. a campaign or day
#'           column) over which the precision is estimated. Default
#'           \code{NULL}: one value for the whole record. Ignored for
#'           \code{method = "allan"}, which is always per \code{UniqueID}.
#' @param tol numerical; for \code{method = "mad"}, the fractional change in
#'            the spacing of \code{POSIX.time} that starts a new run of
#'            constant logging interval (default 0.2, i.e. 20 \%). See Details.
#'
#' @details
#' To use an empirical precision for the minimal detectable flux, pass it to
#' the \code{prec} argument of the import functions (e.g.
#' \code{\link[goFlux]{import.LI7810}}) or write it into the \code{*_prec}
#' columns of the data frame before \code{\link[goFlux]{goFlux}}: the
#' \code{prec} column of the output and \code{\link[goFlux]{MDF}} then use it.
#' See \code{\link[goFlux]{MDF}} for the confidence level (\code{conf}).
#'
#' \strong{Logging interval.} The MAD of first differences describes the noise
#' at the interval the data were logged at, so a record whose interval was
#' changed mid-way (1 s to 10 s, say) would mix two noise scales in one MAD and
#' the result would describe neither. When \code{dataframe} contains
#' \code{POSIX.time}, \code{method = "mad"} therefore splits each record (or
#' \code{by} group) into runs of constant logging interval: a run breaks where
#' the spacing changes by more than \code{tol}; duplicate or backwards
#' timestamps are not counted as intervals; runs shorter than 3 differences are
#' dropped (a lone gap between files is not a run); the differences of all runs
#' with the same interval (rounded to 0.5 s) are pooled into one MAD. The
#' interval is returned in \code{dt_s}. With a single interval the output has
#' one row per record or group, as before. When more than one interval is
#' present, one row \emph{per interval} is returned (so a record or group may
#' occupy several rows) and a warning is given; pick the row whose
#' \code{dt_s} matches the closures you are evaluating. Without
#' \code{POSIX.time} the MAD is taken over the whole record and \code{dt_s} is
#' \code{NA}. To compare an empirical precision with a datasheet value quoted
#' at 1 s, rescale the latter with \code{\link[goFlux]{spec.at.interval}}.
#'
#' @returns A data.frame with one row per record, group or \code{UniqueID}:
#'          the grouping column (if any) and \code{prec} in the units of
#'          \code{gastype}, plus \code{n} (differences used). For
#'          \code{method = "mad"} a trailing column \code{dt_s} gives the
#'          logging interval (seconds) the precision refers to, and a record
#'          with several intervals gets one row per interval (see Details).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{MDF}}, \code{\link[goFlux]{goFlux}},
#'          \code{\link[goFlux]{best.flux}}, \code{\link[goFlux]{spec.at.interval}}
#'
#' @examples
#' data(manID.UGGA)
#' # one robust value for the whole record (dt_s = 1: logged at 1 Hz)
#' empirical.prec(manID.UGGA, "CH4dry_ppb")
#' # per measurement, Allan deviation on the flagged window
#' empirical.prec(manID.UGGA, "CH4dry_ppb", method = "allan")
#' # a record whose logging interval changed from 1 s to 10 s: one row per
#' # interval, with a warning
#' t0 <- as.POSIXct("2024-06-01", tz = "UTC")
#' mixed <- data.frame(POSIX.time = c(t0 + 0:599, t0 + 600 + 10 * (0:299)),
#'                     CH4dry_ppb = 2000 + c(rnorm(600, 0, 1), rnorm(300, 0, 0.4)))
#' empirical.prec(mixed, "CH4dry_ppb")
#' @export
empirical.prec <- function(dataframe, gastype, method = c("mad", "allan"), by = NULL,
                           tol = 0.2) {

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

  x <- dataframe[[gastype]]

  if(method == "mad"){
    has.time <- any(grepl("\\<POSIX.time\\>", names(dataframe)))
    p <- if(has.time) as.numeric(dataframe$POSIX.time) else NULL

    # MAD of first differences per run of constant logging interval.
    # Returns a data.frame with one row per interval: prec, n (differences
    # used) and dt_s (interval, s). Without timestamps: one pooled MAD, dt_s NA.
    f <- function(v, tm){
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
      rownames(out) <- NULL
      out
    }
    if(is.null(by)){
      out <- f(x, p)
      if(nrow(out) > 1){
        warning("record has ", nrow(out), " logging intervals (",
                paste(out$dt_s, collapse = ", "), " s); precision is reported ",
                "per interval (one row each)", call. = FALSE)}
    } else {
      g <- as.character(dataframe[[by]])
      idx.g <- split(seq_along(x), g)
      r <- lapply(names(idx.g), function(k){
        i <- idx.g[[k]]
        o <- f(x[i], if(is.null(p)) NULL else p[i])
        if(nrow(o) > 1){
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
    f <- function(v){ v <- v[!is.na(v)]; dd <- diff(v)
      if(length(dd) < 2) return(c(NA_real_, length(dd)))
      c(stats::sd(dd) / sqrt(2), length(dd)) }
    r <- t(sapply(split(d[[gastype]], as.character(d$UniqueID)), f))
    out <- data.frame(UniqueID = rownames(r), prec = r[, 1], n = as.integer(r[, 2]),
                      row.names = NULL, stringsAsFactors = FALSE)
  }
  out
}
