# Helper functions

# Calculate absolute angle
abs_angle <- function(x, y) {
  dx <- diff(x)
  dy <- diff(y)
  abs <- (pi/2 - atan2(dy, dx)) %% (2*pi)
  return(abs)
}

# Append to HS formula
full_formula <- function(f) {
  reformulate(c(as.character(f[2]), 
                "sl_ + log(sl_) + cos(ta_) + strata(step_id_)"),
              response = "case_")
}

# Function to simulate
sim_track <- function(w, tent, start_loc, start_t,
                      n_avail, T, delta_t = "1 hour") {

  # Create template data.frame
  d <- data.frame(x = rep(NA, T),
                  y = NA,
                  t = seq(from = start_t,
                          by = delta_t,
                          length.out = T))
  # Fill in start location
  d$x[1] <- start_loc[1]
  d$y[1] <- start_loc[2]

  # Generate first step
  # Assume first step goes straight north by mean tentative step length
  d$x[2] <- d$x[1]
  d$y[2] <- d$y[1] + (tent$shp * tent$scl)

  # Generate remaining steps
  for (t in 3:T) {
    # Print status
    cat("\rLocation", t, "of", T, "         ")
    # Absolute angle of previous step
    abs_prev <- abs_angle(x = d$x[(t-2):(t-1)],
                          y = d$y[(t-2):(t-1)])
    # Draw step lengths
    sls <- rgamma(n = n_avail, shape = tent$shp, scale = tent$scl)
    # Draw turn angles
    suppressWarnings({
      ta_circ <- circular::rvonmises(n = n_avail, mu = tent$mu, kappa = tent$kappa)
      tas <- ifelse(ta_circ > pi, ta_circ - 2*pi, ta_circ)
    })

    # Figure out where it lands
    abs_angles <- (abs_prev + tas) %% (2*pi)
    dx <- cos(abs_angles) * sls
    dy <- sin(abs_angles) * sls
    cand_x <- d$x[t-1] + dx
    cand_y <- d$y[t-1] + dy

    # Extract w(x) for each candidate step endpoint
    ws <- extract(w, cbind(cand_x, cand_y))
    # Set any NAs to 0 (hopefully unlikely)
    ws[is.na(ws)] <- 0
    # Normalize
    weight <- ws/sum(ws)

    # Sample step
    step_id <- sample.int(length(weight), size = 1, prob = weight)

    # Keep new location
    d$x[t] <- cand_x[step_id]
    d$y[t] <- cand_y[step_id]
  }
  return(d)
}

# Simulate UD
sim_od1 <- function(object, hs_form, hab, start_loc, start_t, 
                    n_avail, T, delta_t = "1 hour",
                    param_uncert = FALSE) {
  
  if (param_uncert) {
    # Sample parameters
    bb <- MASS::mvrnorm(n = 1, mu = coef(object), 
                        Sigma = stats::vcov(object$model))
  } else {
    bb <- coef(object)
  }
  
  # Indices for habitat selection parms only
  inds <- which(!(grepl("sl_", names(bb))|grepl("ta_", names(bb))))
  # Keep only HS parms
  bb_hs <- bb[inds]
  
  # Tentative movement distributions
  tent <- list(shp = object$sl_$params$shape,
               scl = object$sl_$params$scale,
               mu = object$ta_$params$mu,
               kappa = object$ta_$params$kappa)
  
  # Simulate
  locs <- sim_issf(f = hs_form, b = bb_hs, hab_rast = hab,
                   tent_parm = tent, T = T, n_avail = n_avail)
  
  # Calculate OD
  OD <- locs %>% 
    make_track(x, y, t, crs = 32612) %>% 
    od(trast = hab[[1]])
  
  # Return
  return(OD)
}

sim_od <- function(n_sim, object, hs_form, hab, start_loc, start_t, 
                   n_avail, T, delta_t = "1 hour",
                   param_uncert = TRUE) {
  # Replicate simulation
  ll <- replicate(n = n_sim, expr = {
    sim_od1(object = object, hab = hab, hs_form = hs_form,
            start_loc = start_loc, start_t = start_t,
            n_avail = n_avail, T = T, delta_t = delta_t,
            param_uncert = param_uncert)
  })
  
  # Stack
  ss <- stack(ll)
  
  # Get mean
  mm <- mean(ss)
  
  # Normalize
  MM <- mm/cellStats(mm, sum)
  
  # Return
  return(MM)
}

sim_issf <- function(f, b, hab_rast, tent_parm, T, n_avail) {
  
  ## Pre-compute w(x)
  # Can pre-compute movement-free selection kernel since they are not
  # a function of current location (in this case).
  # Data.frame
  hab_df <- as.data.frame(hab_rast) %>% 
    mutate(cover = factor(cover))
  # Create model matrix
  mm <- model.matrix(f, data = hab_df)
  # Get rid of intercept
  mm <- mm[, 2:ncol(mm), drop = FALSE]
  vv <- exp(mm %*% b)
  # Put in raster
  w <- hab_rast[[1]]
  values(w) <- vv
  names(w) <- "w"
  
  # Simulate
  locs <- sim_track(w = w, tent = tent_parm, 
                    start_loc = c(446000, 4627000), 
                    start_t = ymd_hm("2022-05-01 00:00"), 
                    n_avail = n_avail, T = T, delta_t = "1 hour")
  
  return(locs)
  
}

# Metrics ----
# Concordance
concordance <- function(x) {
  survival::concordance(x$model)$concordance
}

# Coefficient of Determination
cod <- function(x){
  survMisc::rsq(x$model)$cod
}

# Measure of Explained Randomness
mer <- function(x){
  survMisc::rsq(x$model)$mer
}

# Measure of Explained Variation
mev <- function(x){
  survMisc::rsq(x$model)$mev
}

# Scaled rank of used step
srus <- function(x) {
  # Get data
  D <- x$model$model
  # Predict "risk"
  D$risk <- predict(x$model, type = "risk")
  # Easier column to find used
  D$used <- as.character(D[, 1]) == "1"
  # Split data by stratum
  D_sp <- split(D, D$`strata(step_id_)`)
  # Rank
  D_sp <- lapply(D_sp, function(dd) {
    dd$rank <- rank(dd$risk)
    return(dd)
  })
  # Used ranks
  used_ranks <- unlist(
    lapply(D_sp, function(x) {
      x$rank[which(x$used)]
    })
  )
  # Total steps (in case some steps had NAs)
  tot_steps <- unlist(lapply(D_sp, nrow))
  
  # Metric to return
  scaled_rank <- used_ranks/tot_steps
  
  # Return
  return(mean(scaled_rank))
}

rss_v_top <- function(x) {
  # Get data
  D <- x$model$model
  # Predict "risk"
  D$risk <- predict(x$model, type = "risk")
  # Easier column to find used
  D$used <- as.character(D[, 1]) == "1"
  # Split data by stratum
  D_sp <- split(D, D$`strata(step_id_)`)
  # Risk prediction for used 
  risk_used <- unlist(
    lapply(D_sp, function(x) {
      x$risk[which(x$used)]
    })
  )
  # Risk prediction for max
  risk_max <- unlist(
    lapply(D_sp, function(x) {
      # Add [1] in case there are ties for the top risk
      x$risk[which(x$risk == max(x$risk))[1]]
    })
  )
  # RSS
  rss <- risk_used/risk_max
  # Return
  return(mean(rss))
}

# Observed vs Simulated using BA
os_ba <- function(obs, sim) {
  # Normalize
  o <- values(obs)/sum(values(obs))
  s <- values(sim)/sum(values(sim))
  # Calculate overlap (Bhattacharyya's Affinity)
  ba <- amt:::vol_base(o, s, type = "ba")
  # Return
  return(ba)
}

# Observed vs simulated using Spearman's R
os_r <- function(obs, sim) {
  # Calculate Spearman's correlation
  sp_r <- cor(values(obs), values(sim), method = "spearman")
  # Return
  return(sp_r)
}



