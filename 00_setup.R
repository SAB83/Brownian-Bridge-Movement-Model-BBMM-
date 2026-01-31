suppressPackageStartupMessages({
  library(yaml)
  source("R/helpers.R")
})

cfg <- yaml::read_yaml("config/config.yml")

dir_create(
  cfg$outputs$dir,
  file.path(cfg$outputs$dir, "intermediate"),
  file.path(cfg$outputs$dir, "final"),
  file.path(cfg$outputs$dir, "figures")
)

TZ <- cfg$options$timezone %||% "UTC"
CHUNK_SIZE <- cfg$options$chunk_size %||% 5
SEED <- cfg$options$set_seed %||% 1

message("Setup ok. TZ=", TZ, " CHUNK_SIZE=", CHUNK_SIZE, " SEED=", SEED)
