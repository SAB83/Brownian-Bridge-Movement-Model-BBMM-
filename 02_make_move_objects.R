# 02_make_move_objects.R
# Build move objects and reproject to meters for BBMM.
#
# Input: cfg$outputs$cleaned_csv
# Output: cfg$outputs$move_rds

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(sp)
  library(rgdal)
  library(move)
})

dat <- readr::read_csv(cfg$outputs$cleaned_csv, show_col_types = FALSE)
require_cols(dat, c("ID","Long","Lat","LocT"), "cleaned geolocator")

data_list <- split(dat, dat$ID)
message("Tracks: ", length(data_list))

create_move_object <- function(ind_data, proj_out) {
  if (nrow(ind_data) <= 1) return(NULL)

  mv <- move::move(
    x = ind_data$Long,
    y = ind_data$Lat,
    time = as.POSIXct(ind_data$LocT, tz = TZ),
    data = ind_data,
    proj = sp::CRS("+proj=longlat +datum=WGS84"),
    animal = ind_data$ID[1]
  )

  sp::spTransform(mv, sp::CRS(proj_out))
}

proj_out <- cfg$bbmm$projection
move_list <- lapply(data_list, create_move_object, proj_out = proj_out)
move_list <- move_list[!vapply(move_list, is.null, logical(1))]

message("Move objects created: ", length(move_list))
saveRDS(move_list, cfg$outputs$move_rds)
message("Wrote: ", cfg$outputs$move_rds)
