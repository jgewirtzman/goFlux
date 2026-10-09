#' Constant offset between field-notebook closure times and analyzer time
#'
#' Field notebooks record closure start times on a watch or phone; the
#' analyzer stamps its data with its own clock, which may be minutes off (time
#' zone, drift, a reset). \code{find.clock.offset} finds the constant offset
#' \code{delta} such that \code{notebook + delta} lines up with the onsets of
#' concentration rise in the analyzer record.
#'
#' @param dataframe data.frame; a whole imported record (output of an import
#'                  function, e.g. \code{\link[goFlux]{import.UGGA}}) with the
#'                  columns \code{POSIX.time} and \code{gastype}.
#' @param notebook.times POSIXct vector of notebook closure start times, or a
#'                       data.frame (an auxfile) with a POSIXct column
#'                       \code{start.time}.
#' @param gastype character string; the gas used for the alignment. Default
#'                \code{"CO2dry_ppm"} (a sealed chamber on a respiring surface
#'                always makes CO2 rise at closure).
#' @param search numerical vector of length 2; the range of offsets searched
#'               (seconds, analyzer time minus notebook time). Default
#'               \code{c(-1800, 1800)}.
#' @param step numerical; search step (seconds). Default 1.
#' @param window numerical; length (seconds) of the averaging windows before
#'               and after a candidate start. Default 60.
#' @param plot logical; draw the score curve and the aligned traces? Default
#'             \code{TRUE}.
#' @param auxfile data.frame; optional auxfile whose \code{start.time} is
#'                shifted by the best offset (the original is kept as
#'                \code{start.time_notebook}).
#'
#' @details
#' For every candidate offset the function computes, at each notebook start
#' shifted by that offset, an onset score: the rise from the 7 s around that
#' instant to the mean of the \code{window} seconds after it, minus the
#' absolute change from the \code{window} seconds before it. The score is
#' large only where a rise begins from a flat baseline, near zero in the middle
#' of a rise and strongly negative where the chamber is lifted. Each offset is
#' scored by the median over closures, a cross-correlation of the notebook
#' start times with the rise onsets. The whole score curve is returned so that
#' a flat or multi-modal curve can be recognised. The record is interpolated
#' to a 1 s grid first, so any logging interval can be used.
#'
#' @returns A list with \code{offset} (seconds to \emph{add} to the notebook
#'          times), \code{score} (onset score at that offset, units of
#'          \code{gastype}), \code{scores} (data.frame of \code{offset} and
#'          \code{score} over the search range), \code{n.closures} and, when
#'          \code{auxfile} is given, \code{auxfile} with corrected
#'          \code{start.time}.
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{obs.win}}, \code{\link[goFlux]{auto.id.rise}}
#'
#' @examples
#' set.seed(7)
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' tt <- t0 + 0:3599
#' co2 <- 410 + rnorm(3600, 0, 0.5)
#' true.starts <- t0 + c(300, 1200, 2100, 3000)       # analyzer time
#' for (s in c(300, 1200, 2100, 3000))
#'   co2[(s + 1):(s + 300)] <- co2[(s + 1):(s + 300)] + 0.3 * (1:300)
#' record <- data.frame(POSIX.time = tt, CO2dry_ppm = co2)
#' notebook <- true.starts - 137                       # the watch was 137 s behind
#' find.clock.offset(record, notebook, search = c(-600, 600), plot = FALSE)$offset
#' @export
find.clock.offset <- function(dataframe, notebook.times, gastype = "CO2dry_ppm",
                              search = c(-1800, 1800), step = 1, window = 60,
                              plot = TRUE, auxfile = NULL) {

  # Check arguments
  if(!is.data.frame(dataframe)) stop("'dataframe' must be of class data.frame")
  if(!all(c("POSIX.time", gastype) %in% names(dataframe))){
    stop("'dataframe' must contain the columns 'POSIX.time' and '", gastype, "'")}
  nb <- if(is.data.frame(notebook.times)) notebook.times$start.time else notebook.times
  if(!inherits(nb, "POSIXct")){
    stop("'notebook.times' must be POSIXct (or a data.frame with a POSIXct 'start.time')")}
  nb <- nb[!is.na(nb)]
  if(length(nb) == 0) stop("'notebook.times' contains no time")
  if(!is.numeric(search) || length(search) != 2 || search[2] <= search[1]){
    stop("'search' must be two increasing numbers")}
  if(!is.null(auxfile) && !"start.time" %in% names(auxfile)){
    stop("'auxfile' must contain the column 'start.time'")}

  tt <- as.numeric(dataframe$POSIX.time); y <- as.numeric(dataframe[[gastype]])
  ok <- !is.na(tt) & !is.na(y); tt <- tt[ok]; y <- y[ok]
  o <- order(tt); tt <- tt[o]; y <- y[o]
  t0 <- floor(min(tt)); t1 <- ceiling(max(tt))
  grid <- seq(t0, t1, by = 1)
  yg <- stats::approx(tt, y, xout = grid, rule = 1, ties = mean)$y
  has <- !is.na(yg)
  cs <- cumsum(ifelse(has, yg, 0)); cn <- cumsum(has)
  N <- length(grid)
  # mean of yg over grid indices a..b (vectorised, from cumulative sums)
  wmean <- function(a, b) {
    a <- pmax(a, 1); b <- pmin(b, N)
    bad <- a > b
    a2 <- pmin(a, N); b2 <- pmax(b, 1)
    s <- cs[b2] - ifelse(a2 > 1, cs[pmax(a2 - 1, 1)], 0)
    k <- cn[b2] - ifelse(a2 > 1, cn[pmax(a2 - 1, 1)], 0)
    r <- ifelse(k > 0, s / k, NA); r[bad] <- NA; r
  }
  offsets <- seq(search[1], search[2], by = step)
  nbn <- as.numeric(nb)
  score <- vapply(offsets, function(d) {
    idx <- round(nbn + d - t0) + 1
    ctr <- wmean(idx - 3, idx + 3)
    onset <- (wmean(idx, idx + window) - ctr) - abs(ctr - wmean(idx - window, idx))
    if(all(is.na(onset))) NA_real_ else stats::median(onset, na.rm = TRUE)
  }, numeric(1))
  if(all(is.na(score))){
    stop("the notebook times never fall inside the record for any offset in 'search'")}
  best <- offsets[which.max(score)]
  out <- list(offset = best, score = max(score, na.rm = TRUE),
              scores = data.frame(offset = offsets, score = score),
              n.closures = length(nb))
  if(!is.null(auxfile)){
    auxfile$start.time_notebook <- auxfile$start.time
    auxfile$start.time <- auxfile$start.time + best
    out$auxfile <- auxfile
  }
  if(isTRUE(plot)){
    op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1)); on.exit(par(op))
    graphics::plot(offsets, score, type = "l", xlab = "offset added to notebook time (s)",
                   ylab = paste0("median onset score (", gastype, ")"),
                   main = sprintf("best offset = %g s (n = %d closures)", best, length(nb)))
    graphics::abline(v = best, col = "red"); graphics::abline(h = 0, col = "grey60", lty = 2)
    show <- nb[seq_len(min(6, length(nb)))]
    rel <- seq(-2 * window, 4 * window)
    ymat <- sapply(show, function(s) yg[pmin(pmax(round(as.numeric(s) + best - t0) + 1 + rel, 1), N)])
    graphics::matplot(rel, ymat, type = "l", lty = 1, xlab = "seconds from corrected notebook start",
                      ylab = gastype, main = "traces aligned at notebook start + offset")
    graphics::abline(v = 0, col = "red")
  }
  out
}
