# 05_plot_ud_maps.R
# Save individual and composite UD maps to PNG.
#
# Input: cfg$outputs$composite_ud_rds
# Output: PNGs in outputs/figures/

source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(raster) })

obj <- readRDS(cfg$outputs$composite_ud_rds)
ud_resampled <- obj$ud_resampled
composite_ud <- obj$composite_ud

fig_dir <- file.path(cfg$outputs$dir, "figures")
dir_create(fig_dir)

levels <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5)

for (i in seq_along(ud_resampled)) {
  png(file.path(fig_dir, sprintf("ud_individual_%03d.png", i)), width=1400, height=1100, res=150)
  plot(ud_resampled[[i]], main = paste("UD for Track", i))
  contour(ud_resampled[[i]], add=TRUE, levels=levels, col="black")
  dev.off()
}

png(file.path(fig_dir, "ud_composite.png"), width=1400, height=1100, res=150)
plot(composite_ud, main="Composite UD")
contour(composite_ud, add=TRUE, levels=levels, col="black")
dev.off()

message("Saved figures to: ", fig_dir)
