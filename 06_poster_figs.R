# Poster figures

# Load packages ----
library(raster)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ragg)

# Source functions ----
source("99_fun.R")

# Load data ----
# Locations
locs <- read.csv("dat/simulated_data.csv") %>% 
  mutate(t = ymd_hms(t))
# Models
mets <- readRDS("out/metrics.rds")
# Habitat
hab <- stack("dat/habitat_scaled.tif")
names(hab) <- c("forage", "pred", "cover", "dist_to_cent", "wrong")

# Simulated Habitat ----
hab_df <- as.data.frame(hab, xy = TRUE)

# ... forage ----
(forage_plot <- ggplot(hab_df, aes(x = x, y = y, fill = forage)) +
   geom_raster() +
   coord_sf(expand = FALSE) +
   xlab(NULL) +
   ylab(NULL) +
   scale_fill_viridis_c(name = "Forage", option = "E") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank(),
         axis.ticks = element_blank()) +
   NULL)

ggsave("fig/poster/hab_forage.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw")

# ... pred ----
(pred_plot <- ggplot(hab_df, aes(x = x, y = y, fill = pred)) +
   geom_raster() +
   coord_sf(expand = FALSE) +
   xlab(NULL) +
   ylab(NULL) +
   scale_fill_viridis_c(name = "pred", option = "E") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank(),
         axis.ticks = element_blank()) +
   NULL)

ggsave("fig/poster/hab_pred.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw")

# ... cover ----
(cover_plot <- ggplot(hab_df, aes(x = x, y = y, fill = cover)) +
   geom_raster() +
   coord_sf(expand = FALSE) +
   xlab(NULL) +
   ylab(NULL) +
   scale_fill_viridis_c(name = "cover", option = "E") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank(),
         axis.ticks = element_blank()) +
   NULL)

ggsave("fig/poster/hab_cover.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw")

# ... dist_to_cent ----
(dist_to_cent_plot <- ggplot(hab_df, aes(x = x, y = y, fill = dist_to_cent)) +
   geom_raster() +
   coord_sf(expand = FALSE) +
   xlab(NULL) +
   ylab(NULL) +
   scale_fill_viridis_c(name = "dist_to_cent", option = "E") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank(),
         axis.ticks = element_blank()) +
   NULL)

ggsave("fig/poster/hab_dist_to_cent.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw")

# ... wrong ----
(wrong_plot <- ggplot(hab_df, aes(x = x, y = y, fill = wrong)) +
   geom_raster() +
   coord_sf(expand = FALSE) +
   xlab(NULL) +
   ylab(NULL) +
   scale_fill_viridis_c(name = "wrong", option = "E") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank(),
         axis.ticks = element_blank()) +
   NULL)

ggsave("fig/poster/hab_wrong.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw")

# Simulated Tracks ----
bb <- as.vector(extent(hab))

# ... beta1 ----
track_beta1 <- locs %>% 
  filter(sim == "beta1") %>% 
  transmute(x1 = x, y1 = y, t1 = t, 
            x2 = lead(x), y2 = lead(y), t2 = lead(t)) %>% 
  ggplot(aes(x = x1, y = y1)) +
  geom_segment(aes(xend = x2, yend = y2), size = 0.2, alpha = 0.5) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_sf(xlim = bb[1:2], ylim = bb[3:4]) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "green", color = "green")) +
  NULL

ggsave("fig/poster/track_beta1.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw", device = agg_tiff, bg = "green")

# ... beta2 ----
track_beta2 <- locs %>% 
  filter(sim == "beta2") %>% 
  transmute(x1 = x, y1 = y, t1 = t, 
            x2 = lead(x), y2 = lead(y), t2 = lead(t)) %>% 
  ggplot(aes(x = x1, y = y1)) +
  geom_segment(aes(xend = x2, yend = y2), size = 0.2, alpha = 0.5) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_sf(xlim = bb[1:2], ylim = bb[3:4]) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "green", color = "green")) +
  NULL

ggsave("fig/poster/track_beta2.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw", device = agg_tiff, bg = "green")

# ... beta3 ----
track_beta3 <- locs %>% 
  filter(sim == "beta3") %>% 
  transmute(x1 = x, y1 = y, t1 = t, 
            x2 = lead(x), y2 = lead(y), t2 = lead(t)) %>% 
  ggplot(aes(x = x1, y = y1)) +
  geom_segment(aes(xend = x2, yend = y2), size = 0.2, alpha = 0.5) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_sf(xlim = bb[1:2], ylim = bb[3:4]) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "green", color = "green")) +
  NULL

ggsave("fig/poster/track_beta3.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw", device = agg_tiff, bg = "green")

# ... beta4 ----
track_beta4 <- locs %>% 
  filter(sim == "beta4") %>% 
  transmute(x1 = x, y1 = y, t1 = t, 
            x2 = lead(x), y2 = lead(y), t2 = lead(t)) %>% 
  ggplot(aes(x = x1, y = y1)) +
  geom_segment(aes(xend = x2, yend = y2), size = 0.2, alpha = 0.5) +
  geom_point(size = 0.4, alpha = 0.5) +
  coord_sf(xlim = bb[1:2], ylim = bb[3:4]) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "green", color = "green")) +
  NULL

ggsave("fig/poster/track_beta4.tiff", width = 3, height = 2.5, units = "in",
       dpi = 300, compression = "lzw", device = agg_tiff, bg = "green")

# Results ----
res_fig <- mets %>% 
  pivot_longer(concord:os_r,
               names_to = "metric",
               values_to = "fit") %>% 
  mutate(metric = factor(metric, 
                         levels = names(mets)[3:ncol(mets)],
                         labels = c("Concord", "MER", "OOS-C", "RSS v Top", 
                                    "ODBA",
                                    "ODR"))) %>% 
  ggplot(aes(x = beta, y = fit, color = model)) +
  facet_wrap(~ metric) +
  geom_line(size = 1) +
  xlab(expression(beta)) +
  ylab("Metric") +
  scale_color_brewer(name = "Model",
                       breaks = levels(mets$model),
                       labels = c("Full", "Rdcd1", "Rdcd2", 
                                  "Rdcd3", "Wrong"),
                       palette = "BrBG") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 25),
        plot.background = element_rect(fill = "green", color = "green")) +
  NULL

ggsave("fig/poster/main_result.tiff", plot = res_fig,
       width = 10, height = 8, units = "in",
       dpi = 300, compression = "lzw", device = agg_tiff, bg = "green")
