# 01_read_clean_geolocator.R
# Read + clean geolocator fixes and standardize column names for BBMM.
#
# Inputs: cfg$paths$input_csvs
# Output: cfg$outputs$cleaned_csv

source("scripts/00_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(lubridate)
})

input_csvs <- cfg$paths$input_csvs
if (is.null(input_csvs) || length(input_csvs) == 0) stop("No input_csvs in config.")

# Read and stack
dfs <- lapply(input_csvs, function(p) {
  message("Reading: ", p)
  readr::read_csv(p, show_col_types = FALSE) %>% as.data.frame()
})
data_raw <- bind_rows(dfs) %>% standardize_cols()

# Candidate columns (edit if needed)
lat_candidates  <- c("location.lat","location_lat","lat","latitude")
lon_candidates  <- c("location.long","location_long","long","lon","longitude")
time_candidates <- c("d","timestamp","datetime","date_time","loct","LocT","time")
id_candidates   <- c("ID","id","bird_id","tag","individual")
year_candidates <- c("year","Year")

lat_col  <- find_first(lat_candidates,  names(data_raw))
lon_col  <- find_first(lon_candidates,  names(data_raw))
time_col <- find_first(time_candidates, names(data_raw))
id_col   <- find_first(id_candidates,   names(data_raw))
year_col <- find_first(year_candidates, names(data_raw))

if (is.na(lat_col) || is.na(lon_col)) stop("Could not find lat/lon columns. Update candidates in this script.")
if (is.na(time_col)) stop("Could not find a time column. Update candidates in this script.")
if (is.na(id_col)) stop("Could not find an ID column. Update candidates in this script.")

data <- data_raw %>%
  mutate(
    Long = as_num(.data[[lon_col]]),
    Lat  = as_num(.data[[lat_col]]),
    timestamp = parse_timestamp(.data[[time_col]], tz = TZ),
    ID_base = as.character(.data[[id_col]])
  ) %>%
  filter(!is.na(Long) & !is.na(Lat) & !is.na(timestamp))

# Optional: combine ID and year to separate tracks by year
if (!is.na(year_col)) {
  data <- data %>%
    mutate(year = as.character(.data[[year_col]]),
           ID = paste0(ID_base, "_", year))
} else {
  data <- data %>% mutate(ID = ID_base)
}

# Remove duplicates and sort
data <- data %>%
  distinct(ID, timestamp, .keep_all = TRUE) %>%
  arrange(ID, timestamp) %>%
  mutate(LocT = lubridate::force_tz(timestamp, TZ)) %>%
  select(ID, Long, Lat, LocT, everything())

message("Prepared: ", nrow(data), " rows; ", dplyr::n_distinct(data$ID), " tracks")
print(head(data, 10))
print(table(data$ID))

readr::write_csv(data, cfg$outputs$cleaned_csv)
message("Wrote: ", cfg$outputs$cleaned_csv)
