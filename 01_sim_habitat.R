# Simulate data

# Load packages ----
library(NLMR)
library(raster)
library(dplyr)

# Simulate habitat ----
# Forage biomass, predator density, landcover, distance to center

# ... forage ----
# units are g/m2
# reasonable values are 0 - 1000

forage <- nlm_fbm(ncol = 500, nrow = 500, resolution = 50, fract_dim = 0.6,
        user_seed = 1)
forage <- (forage * 1400) - 400
forage[forage < 0] <- 0
# Reduce intermediate forage
forage <- forage^(2/3) * 10

# ... predation risk ----
# units are predators/100 km2
# reasonable values are 0 - 15 (think large predator like wolf)
# want mean to vary through space 
coords <- as.data.frame(forage, xy = TRUE)
# Create mean with sine waves in x and y coords
sin_x <- sin(((coords[, "x"]/16000) * 2 * pi) * 1.5)
sin_y <- sin(((coords[, "y"]/16000) * 2 * pi) * 1.5)
mu_pred <- (sin_x + sin_y + sin_x * sin_y) * 0.4
# Add noise
noise <- nlm_fbm(ncol = 500, nrow = 500, resolution = 50, fract_dim = 1.2,
                 user_seed = 2)
pred_vals <- mu_pred + values(noise)
# Truncate at 0
pred_vals[pred_vals < 0] <- 0
# Stretch values
pred_vals <- pred_vals * 8.5
# Place in raster
pred <- forage
# Reverse values (trying to get slight + corr with forage)
values(pred) <- rev(pred_vals)

# ... cover ----
# 3 categories
set.seed(3)
cover_vals <- nlm_mosaictess(ncol = 500, nrow = 500, resolution = 1,
                             germs = 300)
cover <- forage
values(cover) <- cut(values(cover_vals), c(-0.1, 0.6, 0.85, 1))

# ... distance to center ----
# Units are m
# Included to keep animal from walking off landscape
dist_to_cent <- nlm_distancegradient(ncol = 500, nrow = 500, resolution = 1,
                                     origin = rep(250, 4),
                                     rescale = FALSE)
# Get distances correct
dist_to_cent <- dist_to_cent * 100

# ... wrong ----
# Simulate a variable that does not affect selection
wrong <- nlm_fbm(ncol = 500, nrow = 500, resolution = 50, fract_dim = 0.6,
                 user_seed = 5)

# Create habitat raster ----
hab_df <- as.data.frame(forage, xy = TRUE) %>% 
  rename(forage = layer) %>% 
  # Center coords
  mutate(x = x - mean(x),
         y = y - mean(y)) %>% 
  # Translate to UTM grid 12
  mutate(x = x + 447000,
         y = y + 4626000) %>% 
  # Add other layers
  mutate(pred = values(pred),
         cover = values(cover),
         dist_to_cent = values(dist_to_cent),
         wrong = values(wrong))

# Convert to raster
hab <- rasterFromXYZ(hab_df, res = c(50, 50), crs = 32612)

# Scale
hab_scaled <- scale(hab)
hab_scaled$cover <- as.factor(hab$cover)

# Save ----
# Habitat
writeRaster(hab, "dat/habitat.tif", overwrite = TRUE)
# Habitat scaled
writeRaster(hab_scaled, "dat/habitat_scaled.tif", overwrite = TRUE)
