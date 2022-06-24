# Fit iSSA

# Load packages ----
library(amt)
library(raster)
library(circular)
library(dplyr)
library(tidyr)
library(lubridate)
library(broom)
library(ggplot2)

# Source functions ----
source("99_fun.R")

# Load data ----
locs <- read.csv("dat/simulated_data.csv") %>% 
  mutate(t = ymd_hms(t))
tent <- read.csv("dat/tentative_parms.csv")
hab <- stack("dat/habitat_scaled.tif")
names(hab) <- c("forage", "pred", "cover", "dist_to_cent", "wrong")

# Model formulas ----

f_full <- ~ forage + pred + cover + dist_to_cent

f_red1 <- ~ forage + pred + dist_to_cent

f_red2 <- ~ forage + dist_to_cent

f_red3 <- ~ dist_to_cent

f_wrong <- ~ wrong

# Model combinations ----
dat <- expand.grid(beta = 1:4,
                   model = c("full", paste0("red", 1:3), "wrong")) %>% 
  # Format as tibble for list column
  as_tibble()

# Add formula
dat$f <- lapply(1:nrow(dat), function(r) {
  get(paste0("f_", dat$model[r]))
})

# Format data ----
dat$steps <- lapply(1:nrow(dat), function(r) {
  # Status
  cat(r, "of", nrow(dat), "\n")
  # Get the row
  RR <- dat[r, ]
  # Get the location data
  loc <- locs %>% 
    filter(sim == paste0("beta", RR$beta)) %>% 
    # Format steps
    make_track(x, y, t, crs = 32612) %>% 
    steps()
  # Return
  return(loc)
})

# Fit model ----
dat$issf <- lapply(1:nrow(dat), function(r) {
  # Status
  cat(r, "of", nrow(dat), "\n")
  # Get the row
  RR <- dat[r, ]
  # Get the steps
  stp <- RR$steps[[1]]
  # Get the model formula
  f <- full_formula(dat$f[[r]])
  # Fit the model
  set.seed(123 + r)
  m <- stp %>%  
    random_steps(n_control = 20) %>% 
    filter(step_id_ > 1) %>% 
    extract_covariates(hab) %>% 
    mutate(cover = factor(cover)) %>% 
    fit_issf(f, model = TRUE)
  # Return
  return(m)
})

# Compare to truth ----
true <- data.frame(term = c("forage", "pred", "cover2", "cover3", "dist_to_cent"),
                   base = c(1, -1, 1, -1, -1))
# ... habitat selection ----
dat$hs_true <- lapply(1:nrow(dat), function(i) {
  b <- broom::tidy(dat$issf[[i]]$model) %>% 
    filter(!grepl("sl_", term),
           !grepl("ta_", term)) %>% 
    mutate(lwr = estimate - 1.96 * std.error,
           upr = estimate + 1.96 * std.error) %>% 
    left_join(true, by = "term") %>% 
    mutate(true = base * dat$beta[[i]])
  return(b)
})

# Save ----
# Data
saveRDS(dat, "out/issf_fits.rds")
