#' Import a folder of raw files into one data frame, file by file
#'
#' An alternative to \code{\link[goFlux]{import2RData}} that finds the files
#' of one instrument recursively (with an optional file-name \code{pattern}),
#' skips empty files, calls the instrument's import function on each file with
#' \code{save = FALSE} (nothing is written to disk), reports a file that fails
#' and continues with the others, and returns the merged data.
#'
#' @param path character string; folder containing the raw files.
#' @param instrument character string; the instrument, i.e. the suffix of an
#'                   import function (\code{"UGGA"} for
#'                   \code{\link[goFlux]{import.UGGA}}, \code{"LI7810"}, ...).
#' @param pattern character string; optional regular expression for the file
#'                names, e.g. \code{"_f\\\\d+\\\\.txt$"} for UGGA flux files.
#' @param recursive logical; search subfolders? Default \code{TRUE}.
#' @param min.size numerical; files smaller than this (bytes) are skipped.
#'                 Default 1.
#' @param merge logical; return one data frame (\code{TRUE}, default) or a
#'              named list (\code{FALSE}).
#' @param ... arguments to the import function (\code{date.format},
#'            \code{timezone}, \code{prec}, \code{keep_all}, ...); only the
#'            arguments it accepts are passed on.
#'
#' @details
#' In \code{import2RData}, an error in one file stops the whole import with
#' the message \code{no restart 'muffleError' found} (the error handler calls
#' a restart that does not exist), and the output is always written to
#' \code{getwd()/RData}. \code{import2RData.flat} catches the error of each
#' file instead and returns the data.
#'
#' @returns A data.frame with a \code{source_file} column, ordered by
#'          \code{POSIX.time} (or a list of data.frames with
#'          \code{merge = FALSE}). The files that failed are in the attribute
#'          \code{"failed"}.
#'
#' @include goFlux-package.R
#'
#' @seealso \code{\link[goFlux]{import2RData}}
#'
#' @examples
#' d <- system.file("extdata", "LI7810", package = "goFlux")
#' imp <- import2RData.flat(d, "LI7810", pattern = "\\.data$", timezone = "UTC")
#' nrow(imp); unique(imp$source_file)
#' @export
import2RData.flat <- function(path, instrument, pattern = NULL, recursive = TRUE,
                              min.size = 1, merge = TRUE, ...) {
  fun.name <- paste0("import.", instrument)
  ns <- asNamespace("goFlux")
  if(!exists(fun.name, envir = ns)) stop("goFlux has no import function '", fun.name, "'")
  fun <- get(fun.name, envir = ns)
  files <- list.files(path, pattern = pattern, recursive = recursive, full.names = TRUE)
  files <- files[!dir.exists(files)]
  files <- files[!grepl("^\\._", basename(files))]
  files <- files[file.size(files) >= min.size]
  if(length(files) == 0) stop("no files found in ", path)
  dots <- list(...)
  args <- dots[intersect(names(dots), names(formals(fun)))]
  if("save" %in% names(formals(fun))) args$save <- FALSE
  out <- list(); failed <- character(0)
  for(f in files){
    # call by name: the import functions inspect match.call()[[1]]
    r <- tryCatch(do.call(fun.name, c(list(inputfile = f), args), envir = ns),
                  error = function(e){ message(basename(f), ": ", conditionMessage(e)); NULL })
    if(is.null(r) || nrow(r) == 0){ failed <- c(failed, f); next }
    r <- as.data.frame(r)
    r$source_file <- basename(f)
    out[[basename(f)]] <- r
  }
  message(length(out), " file(s) imported, ", length(failed), " failed")
  if(!isTRUE(merge)){ attr(out, "failed") <- failed; return(out) }
  if(length(out) == 0) stop("no file could be imported")
  res <- do.call(rbind, out); rownames(res) <- NULL
  res <- res[order(res$POSIX.time), ]
  rownames(res) <- NULL
  attr(res, "failed") <- failed
  res
}
