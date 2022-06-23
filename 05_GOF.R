# Goodness-of-fit

# Load packages ----
library(amt)
library(raster)
library(circular)
library(dplyr)
library(lubridate)
library(survival)
library(survMisc) # for survMisc::rsq()

# Source functions ----
source("99_fun.R")

# Load data ----
# Locations
locs <- read.csv("dat/simulated_data.csv") %>% 
  mutate(t = ymd_hms(t))
# Models
dat <- readRDS("out/issf_fits_ods.rds")
# Habitat
hab <- stack("dat/habitat_scaled.tif")
names(hab) <- c("forage", "pred", "cover", "dist_to_cent", "wrong")

# Metrics ----
# Concordance
dat$concord <- sapply(dat$issf, concordance)

# Coefficient of Determination
dat$cod <- sapply(dat$issf, cod)

# Measure of Explained Randomness
dat$mer <- sapply(dat$issf, mer)

# Measure of Explained Variation
dat$mev <- sapply(dat$issf, mev)

# Average rank of used step
dat$srus <- sapply(dat$issf, srus)

# RSS vs Top
dat$rssvt <- sapply(dat$issf, rss_v_top)

# Observed vs Simulated using BA
dat$os_ba <- mapply(os_ba, dat$obs_od, dat$sim_od)

# Observed vs simulated using Spearman's R
dat$os_r <- mapply(os_r, dat$obs_od, dat$sim_od)

# Save ----
mets <- dat %>% 
  select(beta, model, concord:os_r)
saveRDS(mets, "out/metrics.rds")

# Plots ----
ggplot(mets, aes(x = beta, y = concord, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = cod, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = mer, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = mev, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = srus, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = rssvt, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = os_ba, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

ggplot(mets, aes(x = beta, y = os_r, color = model)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()


#########################################################################X
#########################################################################X
# SCRAP ----

# ... Compare simulated to observed UD  ----
# (BA?) -- See Potts et al. 2022

# What to use for observed UD?
# See Schlagel et al 2019. Method implemented in 'amt'.
OD_obs <- dat %>% 
  filter(case_) %>% 
  as_track() %>% 
  od(trast = hab[[1]])

values(OD_obs) <- values(OD_obs)/cellStats(OD_obs, sum)

# For simulation, should simulate same number of timesteps as observed, 
# many times. Simulation should include parameteric uncertainty.
OD_sim <- sim_od(n_sim = 10, object = issf, hab = hab, 
                 start_loc = c(446000, 4627000), 
                 start_t = ymd_hm("2022-05-01 00:00"), 
                 n_avail = 300, T = 2500, delta_t = "1 hour")

{
  par(mfrow = c(1, 2))
  plot(OD_obs, main = "Obs")
  plot(OD_obs, main = "Sim")
}
dev.off()

# Calculate overlap (Bhattacharyya's Affinity)
(ba <- amt:::vol_base(values(OD_obs), values(OD_sim), type = "ba"))

# Calculate Spearman's correlation
(sp_r <- cor(values(OD_obs), values(OD_sim), method = "spearman"))


# NOTE: Normalize risk across steps to use as likelihood ----
# Normalized risk is the probability that a proposed step is
# sampled as the used step
mod_ll <- dat %>% 
  filter(case_) %>% 
  pull(risk_norm) %>% 
  log() %>% 
  sum()

# Null model is that each proposed step is equally likely
#   - Probability that used step is selected is 1/n_steps
#       - Note: this is not exactly true if some proposed steps have NAs
#   - Likelihood across entire model is (1/n_steps)^n_obs
#       - (where n_obs is # observed steps = # strata)
#   - Log-likelihood is n_obs * log(1/n_steps)
n_obs <- dat %>% 
  filter(case_) %>% 
  nrow()
# Equivalent to:
n_obs <- issf$model$nevent

null_ll <- n_obs * log(1/n_steps)

# Cox-Snell method
# See Eqn 1b in Nagelkerke 1991 (Biometrica)
(PR2_CS <- 1 - exp((-2/n_obs) * (mod_ll - null_ll)))

# Compare with 'mer' here (they are identical):
rsq(issf$model)

# Nagelkerke method
# See Eqn 3 in Nagelkerke 1991 (Biometrica)
(PR2_Nag <- PR2_CS/(1 - (exp(2/n_obs * null_ll))))



# ... Numerical evaluation of UHC plots ----
# Bhattacharyya's Affinity
# Kolmogorov-Smirnov test
