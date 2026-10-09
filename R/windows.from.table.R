#' Observation windows from a table of start and end times
#'
#' When the observation windows are known (from field notes, screenshots or a
#' previous processing round), \code{windows.from.table} applies them without
#' clicking. Each window (instrument time) is translated into a per-measurement
#' \code{deadband} (window start minus \code{start.time}) and
#' \code{max.obs.length} (window end minus \code{start.time}), and
#' \code{\link[goFlux]{crop.meas}} is called once with those values from an
#' auxfile. For a uniform dead band or crop, use \code{crop.meas} directly; for
#' an automatic dead band see \code{\link[goFlux]{auto.deadband}}.
#'
#' @param ow.list list of data.frames; output of \code{\link[goFlux]{obs.win}}
#'                (each with \code{UniqueID}, \code{POSIX.time},
#'                \code{start.time} and \code{end.time} or
#'                \code{obs.length}). A single data.frame is also accepted.
#' @param windows data.frame with the columns \code{UniqueID}, \code{start} and
#'                \code{end}: POSIXct, or character strings parsed in
#'                \code{tz}.
#' @param tz character string; time zone used to parse character
#'           \code{start} and \code{end}. Default \code{"UTC"}.
#' @param warn.length numerical; windows with fewer flagged observations give a
#'                    warning. Default 60.
#'
#' @details
#' The output is that of \code{\link[goFlux]{crop.meas}}: \code{start.time},
#' \code{end.time} and \code{obs.length} describe the window, \code{flag} is 1
#' inside it, \code{Etime} counts from its start, and \code{cham.close} and
#' \code{cham.open} keep the recorded closure start and end. A window may start
#' before the recorded \code{start.time} (e.g. after a clock correction):
#' \code{cham.close} still holds the recorded start. Use
#' \code{seal.time = "cham.close"} in \code{\link[goFlux]{qc.flags}}.
#' Measurements of \code{ow.list} without a row in \code{windows} are dropped
#' with a message.
#'
#' @returns A data.frame, all windows bound together (see Details).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{crop.meas}}, \code{\link[goFlux]{click.peak2}},
#'          \code{\link[goFlux]{auto.id.rise}}
#'
#' @examples
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' ow <- list(data.frame(UniqueID = "A", POSIX.time = t0 + 0:400, start.time = t0,
#'                       obs.length = 300, CO2dry_ppm = 400 + 0.3 * (0:400)))
#' win <- data.frame(UniqueID = "A", start = "2024-06-01 10:01:00",
#'                   end = "2024-06-01 10:05:00")
#' manID <- windows.from.table(ow, win)
#' sum(manID$flag); unique(manID$obs.length)
#' @export
windows.from.table <- function(ow.list, windows, tz = "UTC", warn.length = 60) {

  # Check arguments
  if(!is.data.frame(windows)) stop("'windows' must be of class data.frame")
  req <- c("UniqueID", "start", "end")
  if(!all(req %in% names(windows))){
    stop("'windows' must contain the columns ", paste(req, collapse = ", "))}
  if(!inherits(windows$start, "POSIXct")) windows$start <- as.POSIXct(windows$start, tz = tz)
  if(!inherits(windows$end, "POSIXct")) windows$end <- as.POSIXct(windows$end, tz = tz)
  if(any(is.na(windows$start) | is.na(windows$end))) stop("some 'start' or 'end' values could not be parsed")
  if(any(duplicated(windows$UniqueID))) stop("'UniqueID' must be unique in 'windows'")
  if(any(windows$end <= windows$start)) stop("each 'end' must be after its 'start'")

  d <- if(is.data.frame(ow.list)) as.data.frame(ow.list) else do.call(rbind, lapply(ow.list, as.data.frame))
  for(col in c("UniqueID", "POSIX.time", "start.time")){
    if(!col %in% names(d)) stop("'ow.list' must contain the column '", col, "'")}
  if(!"end.time" %in% names(d)){
    if(!"obs.length" %in% names(d)) stop("'ow.list' must contain 'end.time' or 'obs.length'")
    d$end.time <- d$start.time + d$obs.length}
  ids <- unique(as.character(d$UniqueID))
  for(u in setdiff(ids, windows$UniqueID)) message("No window for ", u, "; dropped")
  d <- d[d$UniqueID %in% windows$UniqueID, ]
  if(nrow(d) == 0) return(data.frame())

  # crop.meas() crops forward from start.time (deadband >= 0): a window that
  # starts before the recorded start is cropped from the window start, and the
  # recorded times are put back into cham.close / cham.open afterwards.
  w <- windows[match(d$UniqueID, windows$UniqueID), ]
  rec.close <- d$start.time; rec.open <- d$end.time
  d$start.time <- pmin(d$start.time, w$start)
  d$cham.close <- NULL; d$cham.open <- NULL
  first <- !duplicated(d$UniqueID)
  aux <- data.frame(UniqueID = d$UniqueID[first],
                    deadband = as.numeric(w$start[first] - d$start.time[first], units = "secs"),
                    max.obs.length = as.numeric(w$end[first] - d$start.time[first], units = "secs"),
                    stringsAsFactors = FALSE)
  out <- as.data.frame(crop.meas(d, auxfile = aux, deadband = "aux", max.obs.length = "aux"))
  i <- match(out$UniqueID, d$UniqueID)
  out$cham.close <- rec.close[i]; out$cham.open <- rec.open[i]
  out$flag <- as.numeric(out$flag)
  nf <- tapply(out$flag, out$UniqueID, sum)
  for(u in names(nf)[nf < warn.length]){
    warning(u, ": only ", nf[[u]], " observations in window (< warn.length = ",
            warn.length, ")", call. = FALSE)}
  rownames(out) <- NULL
  out
}
