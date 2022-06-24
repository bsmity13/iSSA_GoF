# Goodness-of-fit

# Load packages ----
library(amt)
library(circular)
library(raster)
library(dplyr)
library(lubridate)

# Source functions ----
source("99_fun.R")

# Load data ----
# Locations
locs <- read.csv("dat/simulated_data.csv") %>% 
  mutate(t = ymd_hms(t))
# Models
dat <- readRDS("out/issf_fits.rds")
# Habitat
hab <- stack("dat/habitat_scaled.tif")
names(hab) <- c("forage", "pred", "cover", "dist_to_cent", "wrong")

# Create observed OD ----
loc_ods <- locs %>% 
  split(~sim) %>% 
  lapply(function(x) {
    cat(unique(x$sim), "\n")
    x %>% 
      make_track(x, y, t, crs = 32612) %>% 
      od(trast = hab[[1]])
  })

# # SAVE
# saveRDS(loc_ods, "out/loc_ods.rds")

# Create simulated ODs ----
# Takes about 10h with n_sim = 25
system.time({
  sim_ods <- lapply(1:nrow(dat), function(r) {
    cat(r, "of", nrow(dat), "\n")
    return(
      sim_od(n_sim = 25, object = dat$issf[[r]], hs_form = dat$f[[r]], hab = hab, 
             start_loc = c(locs$x[1], locs$y[1]), start_t = locs$t[1], 
             n_avail = 500, T = 5000, delta_t = "1 hour", param_uncert = FALSE)
    )
  })
})

# # SAVE
# saveRDS(sim_ods, "out/sim_ods.rds")
# LOAD
sim_ods <- readRDS("out/sim_ods.rds")

# Attach ODs to dat ----
dat$obs_od <- lapply(1:nrow(dat), function(r) {
  return(loc_ods[[dat$beta[r]]])
})

dat$sim_od <- lapply(1:nrow(dat), function(r) {
  return(sim_ods[[r]])
})

# Save ----
saveRDS(dat, "out/issf_fits_ods.rds")
