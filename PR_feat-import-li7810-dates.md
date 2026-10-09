# feat: import.LI7810(dates = NULL) to subset campaign-long .data files before parsing

**Branch:** `feat/import-li7810-dates` (`R/import.LI7810.R` only)

LI-7810 `.data` files grow to 100-200 MB over a campaign (~1.5 million rows
at 1 Hz); `import.LI7810()` reads and time-parses every row even when only a
few field days are needed. This adds `dates = NULL`: when a character vector
of dates (as written in the file's `DATE` column, e.g. `"2022-12-05"`) is
given, the file is read with `readLines()`, the `DATAH` header row is kept
and only the `DATA` rows whose `DATE` field is in `dates` are passed to
`read.delim(text = ...)`; everything after that is unchanged. R only, no
external tools. With `dates = NULL` the function is byte-for-byte the current
one (`all.equal()` on the bundled `LI7810.data` is `TRUE`); a non-matching
`dates` stops with a clear message. Roxygen `@param` added.

Tested on the bundled example (330 rows) and on a 163 MB campaign file (six
field days out of ~1.5 million rows, seconds instead of minutes).
