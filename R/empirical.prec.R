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
#'
#' @details
#' To use an empirical precision for the minimal detectable flux, pass it to
#' the \code{prec} argument of the import functions (e.g.
#' \code{\link[goFlux]{import.LI7810}}) or write it into the \code{*_prec}
#' columns of the data frame before \code{\link[goFlux]{goFlux}}: the
#' \code{prec} column of the output and \code{\link[goFlux]{MDF}} then use it.
#' See \code{\link[goFlux]{MDF}} for the confidence level (\code{conf}).
#'
#' @returns A data.frame with one row per record, group or \code{UniqueID}:
#'          the grouping column (if any) and \code{prec} in the units of
#'          \code{gastype}, plus \code{n} (differences used).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{MDF}}, \code{\link[goFlux]{goFlux}},
#'          \code{\link[goFlux]{best.flux}}
#'
#' @examples
#' data(manID.UGGA)
#' # one robust value for the whole record
#' empirical.prec(manID.UGGA, "CH4dry_ppb")
#' # per measurement, Allan deviation on the flagged window
#' empirical.prec(manID.UGGA, "CH4dry_ppb", method = "allan")
#' @export
empirical.prec <- function(dataframe, gastype, method = c("mad", "allan"), by = NULL) {

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

  x <- dataframe[[gastype]]

  if(method == "mad"){
    f <- function(v){ v <- v[!is.na(v)]; d <- diff(v)
      if(length(d) < 2) return(c(NA_real_, length(d)))
      c(stats::mad(d, constant = 1.4826) / sqrt(2), length(d)) }
    if(is.null(by)){
      r <- f(x)
      out <- data.frame(prec = r[1], n = as.integer(r[2]))
    } else {
      g <- as.character(dataframe[[by]])
      r <- t(sapply(split(x, g), f))
      out <- data.frame(g = rownames(r), prec = r[, 1], n = as.integer(r[, 2]),
                        row.names = NULL, stringsAsFactors = FALSE)
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
