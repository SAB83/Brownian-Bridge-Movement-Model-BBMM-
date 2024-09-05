# Brownian-Bridge-Movement-Model-BBMM-
#Using geolocator data of birds migration and prepare them for BBMM analysis and finding out their migration route and stopover
library(move)

library(raster)

library(tidyverse)

library(moveHMM)

library(dplyr)

library(sp)

library(rgdal)

library(maps)

library(foreach)

library(doParallel)

library(parallel)

library(momentuHMM)

library(geosphere)

library(lubridate)

 

files <- drive_find(pattern = "spring2.main2.csv")

 

# Clean and validate coordinates

data <- data %>%

  mutate(

    location.long = as.numeric(`location-long`),

    location.lat = as.numeric(`location-lat`)

  ) %>%

  filter(!is.na(location.long) & !is.na(location.lat))

# Combine ID and year columns

data <- data %>%

  mutate(ID = paste0(ID, year))

 

# Remove leading single quote from the timestamp column

data$d <- sub("^'", "", data$d)

 

# Convert timestamp to POSIXct

data$timestamp <- as.POSIXct(data$d, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")

 

# Mutate and select required columns

data <- data %>%

  mutate(

    Long = as.numeric(location.long),

    Lat = as.numeric(location.lat),

    LocT = timestamp

  ) %>%

  dplyr::select(ID, Long, Lat, LocT)

 

# Remove rows with missing values in Long, Lat, or LocT

data <- data %>%

  filter(!is.na(Long) & !is.na(Lat) & !is.na(LocT))

 

# Remove duplicates based on ID and timestamp

data <- data[!duplicated(data[c("ID", "LocT")]), ]

 

# Order data by timestamp

data <- data[order(data$LocT), ]

 

# Ensure all timestamps are in UTC

data$LocT <- force_tz(data$LocT, "UTC")

 

# Verify consistency of timezones

print(unique(sapply(data$LocT, tz)))  # Should print only "UTC"

 

# Check the final prepared data

print(str(data))

print(head(data))

# Further diagnostics: Check unique IDs and counts

print(unique(data$ID))

print(table(data$ID))

# Split data by ID

data_list <- split(data, data$ID)

print(length(data_list))  # Should print the number of unique IDs

print(head(data_list, 5))  # Inspect the first few entries

 

# Function to create a move object and transform it

create_move_object <- function(ind_data) {

  if (nrow(ind_data) > 1) {

    move_obj <- move(x = ind_data$Long, y = ind_data$Lat, time = as.POSIXct(ind_data$LocT),

                     data = ind_data, proj = CRS("+proj=longlat +datum=WGS84"), animal = ind_data$ID[1])

    move_obj <- spTransform(move_obj, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"))

    return(move_obj)

  } else {

    return(NULL)

  }

}

# Apply the function to create move objects

move_list <- lapply(data_list, create_move_object)

# Filter out NULL elements from move_list

move_list <- move_list[!sapply(move_list, is.null)]

print(length(move_list))  # Should print the number of successfully created move objects

 

# Debugging: Print the first few move objects

print(head(move_list, 5))

# Define the parameters

margin.size <- 5

window.size <- 31

location.error <- 8 # in meters

grid.size <- 50 # Increasing to 20 km grid size in meters

 

# Function to create a dBBMM object with increased extent

create_dbbmm <- function(move_obj) {

  if (!is.null(move_obj)) {

    cat("Creating dBBMM for ID:", move_obj@idData$individual.local.identifier[1], "\n")

    tryCatch({

      dbbmm_obj <- brownian.bridge.dyn(object = move_obj, ext = 0.5, location.error = location.error, window.size = window.size)

      return(dbbmm_obj)

    }, error = function(e) {

      cat("Error creating dBBMM for ID:", move_obj@idData$individual.local.identifier[1], "\n")

      print(e)

      return(NULL)

    })

  } else {

    return(NULL)

  }

}

# Apply the function to create dBBMM objects

chunk_size <- 5

ud_list <- list()

move_list <- list()

for (i in seq(1, length(data_list), by = chunk_size)) {

  chunk <- data_list[i:min(i + chunk_size - 1, length(data_list))]

 

  # Create move objects

  move_chunk <- lapply(chunk, create_move_object)

  move_chunk <- move_chunk[!sapply(move_chunk, is.null)]

 

  # Create dBBMM objects

  dbbmm_chunk <- lapply(move_chunk, create_dbbmm)

  dbbmm_chunk <- dbbmm_chunk[!sapply(dbbmm_chunk, is.null)]

 

  # Calculate the UDs

  for (dbbmm_obj in dbbmm_chunk) {

    ud_ind <- getVolumeUD(dbbmm_obj)

    ud_ind <- ud_ind / sum(values(ud_ind), na.rm = TRUE) # Normalize UD

    ud_list[[length(ud_list) + 1]] <- ud_ind

  }

 

  move_list <- c(move_list, move_chunk)

}

# Check the number of successfully created UDs

print(length(ud_list))

# Merge all UDs to create a combined extent

combined_extent <- extent(do.call(merge, lapply(ud_list, extent)))

reference_raster <- raster(

  xmn = xmin(combined_extent),

  xmx = xmax(combined_extent),

  ymn = ymin(combined_extent),

  ymx = ymax(combined_extent),

  nrows = 100,

  ncols = 100,

  crs = CRS("+proj=longlat +datum=WGS84")

)

# Assuming `ud_list` and `reference_raster` are already created as per the previous steps

 

# Resample UDs to the reference raster

ud_list_resampled <- lapply(ud_list, function(ud) {

  resample(ud, reference_raster, method = "bilinear")

})

 

# Sum individual UDs

composite_ud <- Reduce('+', ud_list_resampled)

composite_ud <- composite_ud / sum(values(composite_ud), na.rm = TRUE)  # Normalize

 

# Set up a larger plot window

par(mfrow = c(1, 1), oma = c(0, 0, 2, 0), mar = c(4, 4, 2, 2))  # Set up the plot layout

 

# Define plot limits based on the extent of the reference raster

xlim <- c(min(reference_raster@extent@xmin, na.rm = TRUE), max(reference_raster@extent@xmax, na.rm = TRUE))

ylim <- c(min(reference_raster@extent@ymin, na.rm = TRUE), max(reference_raster@extent@ymax, na.rm = TRUE))

 

# Plot each UD with larger map size and better view

for (i in seq_along(ud_list_resampled)) {

  plot(ud_list_resampled[[i]], main = paste("UD for Individual", i), xlim = xlim, ylim = ylim, col = terrain.colors(100))

  contour(ud_list_resampled[[i]], add = TRUE, levels = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5), col = "black")

}

plot(composite_ud, main = "Composite UD", xlim = xlim, ylim = ylim, col = terrain.colors(100))

contour(composite_ud, add = TRUE, levels = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5), col = "black")

 

# Convert columns to numeric and filter out NA values

data <- data %>%

  mutate(

    location_long = as.numeric( `location-long`),

    location_lat = as.numeric(`location-lat`)

  ) %>%

  filter(!is.na(location_long) & !is.na(location_lat))

 

# Create unique ID by combining ID and year

data <- data %>%

  mutate(ID = paste0(ID, year))

 

# Clean up date column and convert to POSIXct

data$d <- sub("^'", "", data$d)

data$timestamp <- as.POSIXct(data$d, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")

 

# Add new columns and select required columns

data <- data %>%

  mutate(

    Long = as.numeric(location_long),

    Lat = as.numeric(location_lat),

    year = as.factor(year),

    LocT = timestamp

  ) %>%

  dplyr::select(ID, Long, Lat, LocT, year)

 

# Filter and arrange data

data <- data %>%

  filter(!is.na(Long) & !is.na(Lat) & !is.na(LocT)) %>%

  arrange(ID, LocT)

data <- data[!duplicated(data[c("ID", "LocT")]), ]

 

# Ensure all timestamps are in UTC

data$LocT <- force_tz(data$LocT, "UTC")

 

# Split data by ID

data_list <- split(data, data$ID)

print(length(data_list))  # Should print the number of unique IDs

print(head(data_list, 5))  # Inspect the first few entries

 

 

# Function to create a move object and transform it

create_move_object <- function(ind_data) {

  if (nrow(ind_data) > 1) {

    move_obj <- move(x = ind_data$Long, y = ind_data$Lat, time = as.POSIXct(ind_data$LocT),

                     data = ind_data, proj = CRS("+proj=longlat +datum=WGS84"), animal = ind_data$ID[1])

    move_obj <- spTransform(move_obj, CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"))

    return(move_obj)

  } else {

    return(NULL)

  }

}

 
