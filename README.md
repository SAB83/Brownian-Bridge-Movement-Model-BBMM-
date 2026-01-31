# Geolocator dBBMM Migration Pipeline (R)

A clean, GitHub-ready pipeline to prepare bird geolocator tracks and estimate migration routes / utilization distributions (UDs)
using a **dynamic Brownian Bridge Movement Model** (dBBMM) via the `move` package.

**No data are included** (private/non-public geolocator data). You point the pipeline to your local files using `config/config.yml`
(which is git-ignored).

## What the pipeline does (scripts)

1. **01_read_clean_geolocator.R**
   - Reads one or more CSVs
   - Maps your lat/lon/time columns into standard names
   - Converts coordinates to numeric and drops invalid rows
   - Parses timestamps and enforces a single timezone (UTC by default)
   - Creates a per-track identifier (`ID_year` if `year` exists)
   - Removes duplicate fixes by `ID + time`
   - Output: `outputs/intermediate/geolocator_cleaned.csv`

2. **02_make_move_objects.R**
   - Splits cleaned data by `ID`
   - Builds `move` objects (WGS84)
   - Transforms to a **projected CRS in meters** for BBMM
   - Output: `outputs/intermediate/move_objects.rds`

3. **03_fit_dbbmm.R**
   - Fits dBBMM per individual with `move::brownian.bridge.dyn()`
   - Robust error handling + chunk processing
   - Converts each model to a UD raster (`getVolumeUD`) and normalizes
   - Output: `outputs/intermediate/dbbmm_objects.rds`, `outputs/intermediate/ud_list.rds`

4. **04_build_composite_ud.R**
   - Creates a reference raster covering all individuals
   - Resamples all UDs to the same grid
   - Sums + normalizes to create a composite UD
   - Output: `outputs/final/composite_ud.rds`

5. **05_plot_ud_maps.R**
   - Saves PNGs for each individual UD and the composite UD (with contours)
   - Output: `outputs/figures/*.png`

## Install packages

```r
install.packages(c(
  "move","raster","sp","rgdal",
  "dplyr","stringr","readr","lubridate","yaml"
))
```

> If `rgdal` is difficult to install on your system, you can still run most of the pipeline;
the only place it's needed is `spTransform()` in script 02. We can switch that step to `sf` if needed.

## Configure paths

```bash
cp config/config_example.yml config/config.yml
```

Edit `config/config.yml` to point to your local CSV file(s).

## Run

```r
source("scripts/01_read_clean_geolocator.R")
source("scripts/02_make_move_objects.R")
source("scripts/03_fit_dbbmm.R")
source("scripts/04_build_composite_ud.R")
source("scripts/05_plot_ud_maps.R")
```

## Expected input columns

Your CSV should contain (names can vary; script 01 searches common variants):

- ID column (e.g., `ID`)
- Year column (optional; e.g., `year`)
- Latitude column (e.g., `location-lat`)
- Longitude column (e.g., `location-long`)
- Datetime column (e.g., `d`)

## Suggested GitHub repo name

**geolocator-dbbmm-migration-pipeline**
