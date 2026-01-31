# 04_build_composite_ud.R
# Resample UDs to a shared grid and build a composite UD.
#
# Input: cfg$outputs$ud_rds
# Output: cfg$outputs$composite_ud_rds

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(raster) })

ud_list <- readRDS(cfg$outputs$ud_rds)
if (length(ud_list) == 0) stop("ud_list is empty. Run script 03 first.")

combined_extent <- raster::extent(do.call(merge, lapply(ud_list, raster::extent)))

reference_raster <- raster::raster(
  xmn = raster::xmin(combined_extent),
  xmx = raster::xmax(combined_extent),
  ymn = raster::ymin(combined_extent),
  ymx = raster::ymax(combined_extent),
  nrows = cfg$raster$nrows,
  ncols = cfg$raster$ncols,
  crs = raster::crs(ud_list[[1]])
)

method <- cfg$raster$resample_method %||% "bilinear"

ud_resampled <- lapply(ud_list, function(ud) {
  raster::resample(ud, reference_raster, method = method)
})

composite_ud <- Reduce(`+`, ud_resampled)
composite_ud <- composite_ud / sum(raster::values(composite_ud), na.rm = TRUE)

saveRDS(list(reference_raster = reference_raster,
             ud_resampled = ud_resampled,
             composite_ud = composite_ud),
        cfg$outputs$composite_ud_rds)

message("Wrote: ", cfg$outputs$composite_ud_rds)
