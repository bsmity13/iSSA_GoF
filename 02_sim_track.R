# Simulate data

# Load packages ----
library(raster)
library(circular)
library(dplyr)
library(lubridate)

# Source functions ----
source("99_fun.R")

# Load habitat ----
hab <- stack("dat/habitat_scaled.tif")
names(hab) <- c("forage", "pred", "cover", "dist_to_cent", "wrong")

# Create betas ----
# Betas for forage, pred, cover == 2, cover == 3, dist_to_cent, wrong
# Forage always +
# Pred always -
# Cover: 2 > 1 > 3
# Dist_to_cent always -
# Wrong always 0
beta1 <- c(1, -1, 1, -1, -1, 0)
beta2 <- beta1 * 2
beta3 <- beta1 * 3
beta4 <- beta1 * 4

beta <- list(beta1, beta2, beta3, beta4)

# Create selection-free movement distributions ----
tent <- list()

# ... step lengths ----
# Gamma distribution
tent$shp <- 5
tent$scl <- 150
sl <- 0:5000
plot(sl, dgamma(sl, shape = tent$shp, scale = tent$scl), type = "l")

# ... turn angles ----
# von Mises distribution
tent$mu <- 0
tent$kappa <- 1
ta <- seq(-pi, pi, length.out = 200)
plot(ta, dvonmises(ta, mu = tent$mu, kappa = tent$kappa), type = "l",
     ylim = c(0, 0.5))
# How much more likely to go straight than turn around?
suppressWarnings({
  dvonmises(0, mu = tent$mu, kappa = tent$kappa)/
    dvonmises(pi, mu = tent$mu, kappa = tent$kappa)
  })

# Simulate ----
# For each set of betas, simulate dataset under an iSSA

# Number of timesteps
T <- 5000
# Number of candidate steps
n_avail <- 500

# Go!
set.seed(123456)
locs <- lapply(beta, function(B) {
  sim_issf(f = ~ forage + pred + cover + dist_to_cent + wrong, b = B, 
           hab_rast = hab, tent_parm = tent, T = T, n_avail = n_avail)
})

names(locs) <- paste0("beta", 1:4)

# Quick check:
for (i in 1:length(locs)) {
  plot(hab[[1]], main = i)
  points(locs[[i]]$x, locs[[i]]$y, pch = ".")
  lines(locs[[i]]$x, locs[[i]]$y)
}

# Combine
sim_dat <- bind_rows(locs, .id = "sim")

# Save
write.csv(sim_dat, "dat/simulated_data.csv", row.names = FALSE)
write.csv(as.data.frame(tent), "dat/tentative_parms.csv", row.names = FALSE)
