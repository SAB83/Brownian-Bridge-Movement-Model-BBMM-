# 03_fit_dbbmm.R
# Fit dBBMM per track and produce UDs.
#
# Input: cfg$outputs$move_rds
# Outputs: cfg$outputs$dbbmm_rds, cfg$outputs$ud_rds

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(move)
  library(raster)
})

set.seed(SEED)

move_list <- readRDS(cfg$outputs$move_rds)
if (length(move_list) == 0) stop("move_list is empty. Run script 02 first.")

location_error <- cfg$bbmm$location_error_m
window_size <- cfg$bbmm$window_size
ext <- cfg$bbmm$ext

create_dbbmm <- function(move_obj) {
  id <- move_obj@idData$individual.local.identifier[1]
  cat("Creating dBBMM for ID:", id, "\n")
  tryCatch({
    move::brownian.bridge.dyn(
      object = move_obj,
      ext = ext,
      location.error = location_error,
      window.size = window_size
    )
  }, error = function(e) {
    cat("Error for ID:", id, "\n")
    message(e$message)
    NULL
  })
}

ids <- seq_along(move_list)
chunk_size <- CHUNK_SIZE

dbbmm_list <- list()
ud_list <- list()

for (i in seq(1, length(ids), by = chunk_size)) {
  idx <- ids[i:min(i + chunk_size - 1, length(ids))]
  chunk <- move_list[idx]

  dbbmm_chunk <- lapply(chunk, create_dbbmm)
  dbbmm_chunk <- dbbmm_chunk[!vapply(dbbmm_chunk, is.null, logical(1))]

  for (obj in dbbmm_chunk) {
    ud <- move::getVolumeUD(obj)
    ud <- ud / sum(raster::values(ud), na.rm = TRUE)  # normalize
    ud_list[[length(ud_list) + 1]] <- ud
  }

  dbbmm_list <- c(dbbmm_list, dbbmm_chunk)
}

message("dBBMM fits: ", length(dbbmm_list))
message("UD rasters: ", length(ud_list))

saveRDS(dbbmm_list, cfg$outputs$dbbmm_rds)
saveRDS(ud_list, cfg$outputs$ud_rds)
message("Wrote: ", cfg$outputs$dbbmm_rds)
message("Wrote: ", cfg$outputs$ud_rds)
