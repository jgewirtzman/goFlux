#' Fluxes, detection class and quality flags for one gas in one call
#'
#' Runs \code{\link[goFlux]{goFlux}} -> \code{\link[goFlux]{best.flux}} ->
#' \code{\link[goFlux]{flux.class}} -> \code{\link[goFlux]{qc.flags}} (and
#' \code{\link[goFlux]{co2.tracer}} when a CO2 result is given) on the flagged
#' data of all measurements, and records every option used in a
#' \code{settings} list for the methods section. Columns are only appended:
#' the \code{best.flux} columns are returned as \code{best.flux} produced them.
#' Ebullition is not handled here: incubations over water go to
#' \code{\link[goFlux]{goAquaFlux}}.
#'
#' @param dataframe data.frame; flagged concentration data of all
#'                  measurements (output of \code{\link[goFlux]{click.peak2}},
#'                  \code{\link[goFlux]{crop.meas}},
#'                  \code{\link[goFlux]{windows.from.table}} or
#'                  \code{\link[goFlux]{auto.id.rise}}), including the rows
#'                  around the window (\code{flag == 0}) for the ambient check.
#' @param gastype character string; the gas to process (one per call).
#' @param auxfile data.frame; optional, one row per \code{UniqueID}; the
#'                \code{by} column is taken from it when the flux table does
#'                not have it.
#' @param by character string; optional column defining the groups (e.g.
#'           analyzer x field day) of the empirical precision
#'           (\code{\link[goFlux]{flux.class}}) and of the \code{qc.c0}
#'           median (\code{\link[goFlux]{qc.flags}}).
#' @param prec,conf passed to \code{\link[goFlux]{flux.class}}: \code{prec =
#'             NULL} (default) uses the group second-difference precision of
#'             \code{\link[goFlux]{empirical.prec}}; \code{conf = 0.95} gives
#'             z = 1.96 (a benchmark multiplier, see
#'             \code{\link[goFlux]{MDF}}).
#' @param qc named list of arguments to \code{\link[goFlux]{qc.flags}} (e.g.
#'           \code{list(min.obs = 120, leak.rate = 1e-4)}), or \code{FALSE}
#'           to skip the flags. Default \code{list()}: the defaults of
#'           \code{qc.flags}.
#' @param co2.flux.result data.frame; optional \code{best.flux} output for
#'                        \code{"CO2dry_ppm"} of the same measurements, for
#'                        the \code{co2.tracer} column. Do not use on open
#'                        water, lit foliage or dead wood.
#' @param best.flux.args named list of arguments to
#'                       \code{\link[goFlux]{best.flux}}.
#' @param ... further arguments to \code{\link[goFlux]{goFlux}} (e.g.
#'            \code{H2O_col}, \code{k.min}, \code{prec}).
#'
#' @details
#' The precision given to \code{goFlux} (its \code{prec} argument or the
#' \code{*_prec} columns) sets goFlux's own \code{MDF} column and the kappa-max
#' bound of the HM model, and so can change the fluxes; the precision used by
#' \code{flux.class} (\code{det.prec}) only sets \code{det.MDF} and
#' \code{det.class}. To make them the same, compute the precision first
#' (\code{\link[goFlux]{empirical.prec}}), write it into the \code{*_prec}
#' columns, and pass it as \code{prec} too.
#'
#' @returns A list with \code{fluxes} (one row per \code{UniqueID}),
#'          \code{dataframe} (the data fitted), \code{gastype} and
#'          \code{settings} (every option used, the goFlux version and, when
#'          installed from GitHub, the commit).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{write.outputs}}
#'
#' @examples
#' data(manID.UGGA)
#' res <- process.fluxes(manID.UGGA, "CH4dry_ppb")
#' res$fluxes[, c("UniqueID", "best.flux", "MDF", "det.MDF", "det.class", "qc.any")]
#' str(res$settings, max.level = 1)
#' @export
process.fluxes <- function(dataframe, gastype, auxfile = NULL, by = NULL,
                           prec = NULL, conf = 0.95, qc = list(),
                           co2.flux.result = NULL, best.flux.args = list(), ...) {

  # Check arguments
  if(!is.data.frame(dataframe)) stop("'dataframe' must be of class data.frame")
  if(!gastype %in% names(dataframe)) stop("'dataframe' must contain the column '", gastype, "'")
  if(!"flag" %in% names(dataframe)) stop("'dataframe' must contain 'flag' (select the windows first)")
  if(!isFALSE(qc) && !is.list(qc)) stop("'qc' must be a list of qc.flags() arguments or FALSE")
  if(!is.list(best.flux.args)) stop("'best.flux.args' must be a list")

  d <- as.data.frame(dataframe)
  flux <- goFlux(d, gastype, ...)
  best <- as.data.frame(do.call(best.flux, c(list(flux.result = flux), best.flux.args)))
  if(!is.null(by) && !by %in% names(best)){
    src <- if(!is.null(auxfile) && by %in% names(auxfile)) auxfile else d
    if(!by %in% names(src)) stop("column '", by, "' not found in 'auxfile' or 'dataframe'")
    best[[by]] <- src[[by]][match(best$UniqueID, src$UniqueID)]
  }
  best <- flux.class(best, dataframe = d, gastype = gastype, prec = prec, by = by, conf = conf)
  if(!isFALSE(qc)){
    best <- do.call(qc.flags, c(list(flux.result = best, dataframe = d, gastype = gastype,
                                     by = by), qc))
  }
  if(!is.null(co2.flux.result)) best$co2.tracer <- co2.tracer(co2.flux.result, best)

  desc <- utils::packageDescription("goFlux")
  settings <- list(
    goFlux.version = as.character(utils::packageVersion("goFlux")),
    goFlux.commit = if(!is.null(desc$RemoteSha)) desc$RemoteSha else NA_character_,
    gastype = gastype, by = by,
    prec = if(is.null(prec)) "empirical.prec(method = 'hadamard') per group" else
      if(is.data.frame(prec)) "supplied per UniqueID" else prec,
    conf = conf,
    det.prec.median = stats::median(best$det.prec, na.rm = TRUE),
    qc = if(isFALSE(qc)) NULL else qc,
    co2.tracer = !is.null(co2.flux.result),
    goFlux.args = list(...), best.flux.args = best.flux.args,
    n.measurements = nrow(best),
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  list(fluxes = best, dataframe = d, gastype = gastype, settings = settings)
}

#' Write the results of process.fluxes() to a folder
#'
#' Writes \code{fluxes.csv} (one row per measurement), \code{settings.json}
#' (every option used, for the methods section) and, with \code{plots = TRUE},
#' \code{flux_plots.pdf} (\code{\link[goFlux]{flux.plot}} and
#' \code{\link[goFlux]{flux2pdf}}).
#'
#' @param result list; output of \code{\link[goFlux]{process.fluxes}}.
#' @param dir character string; output folder (created if needed).
#' @param plots logical; write the plots? Default \code{TRUE}.
#' @param ... further arguments to \code{\link[goFlux]{flux.plot}}.
#'
#' @returns The paths written, invisibly (named character vector).
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{process.fluxes}}
#'
#' @examples
#' data(manID.UGGA)
#' res <- process.fluxes(manID.UGGA, "CH4dry_ppb")
#' out <- write.outputs(res, file.path(tempdir(), "goFlux_out"), plots = FALSE)
#' basename(out)
#' @export
write.outputs <- function(result, dir, plots = TRUE, ...) {
  if(!is.list(result) || !all(c("fluxes", "dataframe", "gastype", "settings") %in% names(result))){
    stop("'result' must be the output of process.fluxes()")}
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fmt <- function(df){
    df[] <- lapply(df, function(v) if(inherits(v, "POSIXct")) format(v, "%Y-%m-%d %H:%M:%S") else v)
    df }
  paths <- character(0)
  f <- file.path(dir, "fluxes.csv")
  utils::write.csv(fmt(result$fluxes), f, row.names = FALSE); paths["fluxes"] <- f
  f <- file.path(dir, "settings.json")
  jsonlite::write_json(result$settings, f, auto_unbox = TRUE, pretty = TRUE, null = "null",
                       na = "null", digits = NA)
  paths["settings"] <- f
  if(isTRUE(plots)){
    pl <- flux.plot(result$fluxes, result$dataframe, result$gastype, ...)
    f <- file.path(dir, "flux_plots.pdf")
    flux2pdf(pl, outfile = f)
    paths["plots"] <- f
  }
  invisible(paths)
}
