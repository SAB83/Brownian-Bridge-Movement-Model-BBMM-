suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

dir_create <- function(...) {
  for (p in c(...)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

require_cols <- function(df, cols, name="data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop("Missing columns in ", name, ": ", paste(miss, collapse=", "))
  invisible(TRUE)
}

# Convert time strings to POSIXct in a given time zone (default UTC).
# Also removes a leading single quote sometimes produced by spreadsheet exports.
parse_timestamp <- function(x, tz="UTC") {
  x <- as.character(x)
  x <- str_replace(x, "^'+", "")
  x <- str_trim(x)

  ts <- suppressWarnings(as.POSIXct(x, format="%Y-%m-%d %H:%M:%OS", tz=tz))
  if (all(is.na(ts))) ts <- suppressWarnings(as.POSIXct(x, format="%Y/%m/%d %H:%M:%OS", tz=tz))
  if (all(is.na(ts))) ts <- suppressWarnings(as.POSIXct(x, tz=tz))  # last resort
  ts
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

standardize_cols <- function(df) {
  names(df) <- make.names(names(df))
  df
}

find_first <- function(cands, nms) {
  cands <- intersect(cands, nms)
  if (length(cands) == 0) NA_character_ else cands[1]
}
