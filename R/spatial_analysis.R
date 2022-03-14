# SET UP
# Install and load packages
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library e.g. readOGR(), spTransform()
library(sf) # encoding spatial vector data
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # aligning  graphs
library(ggrepel) # for adding labels like geom_text_repel
library(brms) # Bayesian Regression Models
library(rstan)
library(dplyr)
library(tidyr)
library(rnaturalearth) # load map
library(rnaturalearthdata) # load data for map

# mytheme and mycol from species_analysis.R

iucn       <- rbind(c("#D6302D","#D1633A","#CA9D00","#6E6E6E"))
turtle_loc <- clean_meta_between %>%
  subset(!is.na(lon) & !is.na(lat)) %>%
  dplyr::mutate(red_list = factor(red_list, levels = c("Critically endangered", "Endangered", "Vulnerable","Data deficient")))

## PROTECTED BEACHES PER COUNTRY ##-----------------------------------------------
# Download the data
download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", destfile = 'coastlines.zip')

# Load datasets - coastline, country maps, and turtle by country data
coastlines      <- rgdal::readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
coastlines_df   <- sp::SpatialLinesDataFrame(coastlines, coastlines@data) # turn the data into a spatial data frame
world           <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world_turtle_df <- read.csv("world_turtle.csv")
world_merged    <- merge(world, world_turtle_df, by.x = "name", by.y = "name")

# Plot Fig 2a
spatial_plot <- ggplot(data = world_merged) + mytheme +
  geom_sf(aes(fill = all_protected_perc), colour = NA, size = 0.2, alpha = 0.8) + 
  geom_path(data = coastlines_df, aes(x = long, y = lat, group = group), size = 0.1) +
  geom_point(data = turtle_loc, aes(x = lon, y = lat), colour = "black", size = 3) +
  geom_point(data = turtle_loc, aes(x = lon, y = lat, colour = species), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(limits = c(-60, 80), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_colour_manual(values = mycol) +
  colorspace::scale_fill_continuous_sequential(palette = "BluGrn", na.value = "#F1EEE7") +
  theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black"), legend.position = "bottom")

# ANALYSIS  
# Load dataset
worldgroup_dat   <- read.csv("world_turtle_grouped.csv")
worldgroup_dat_2 <- worldgroup_dat %>%
  dplyr::filter(!species %in% c("Lepidochelys_kempii", "Natator_depressus")) %>%
  droplevels()

# GENERAL STAN SPECS
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use

spatial_priors <- c(prior(normal(0, 10), "b"), 
                    prior(normal(0, 10), "Intercept"),
                    prior(student_t(3, 0, 20), "sd"), 
                    prior(student_t(3, 0, 20), "sigma"))

set.seed(10)
protect_mass_model <- brms::brm(mean_mass ~ protected_perc + mean_temp + (1 + protected_perc | species), 
                             data    = worldgroup_dat_2, 
                             family  = gaussian(),
                             prior   = spatial_priors,
                             chains  = 4, 
                             cores   = 4,
                             iter    = 5e3, 
                             warmup  = 2.5e3, 
                             control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(protect_mass_model) # show summary of results

# check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
worldgroup_dat2       <- worldgroup_dat_2 %>% tidyr::drop_na(mean_mass, protected_perc, mean_temp)
protect_y             <- worldgroup_dat2$mean_mass
protect_yrep          <- posterior_predict(protect_mass_model, draws = 500)
protect_species       <- worldgroup_dat2$species
bayesplot::ppc_scatter_avg(protect_y, protect_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + geom_point(aes(colour = protect_species), size = 2.5) + scale_colour_manual(values = mycol)


# Extract fixed and random effect from model
re_model_only <- tidyr::crossing(protected_perc = seq(min(worldgroup_dat_2$protected_perc[!is.na(worldgroup_dat_2$protected_perc)]), 
                                                      max(worldgroup_dat_2$protected_perc[!is.na(worldgroup_dat_2$protected_perc)]),
                                                      length.out = 100),
                                 mean_temp      = mean(worldgroup_dat_2$mean_temp, na.rm = TRUE),
                                 species        = unique(worldgroup_dat_2$species)) %>%
  tidybayes::add_fitted_draws(protect_mass_model,
                              scale = "response", n = 100)

re_model_summary <- re_model_only %>%
  group_by(species, protected_perc) %>%
  summarize(.value = mean(.value))

# Plot Fig 2b
protect_plot <- ggplot() +
  geom_line(data = re_model_only, aes(x = protected_perc, y = .value, group = paste(species, .draw), colour = species), alpha = .1) +
  geom_line(data = re_model_summary, aes(x = protected_perc, y = .value, group = species), size = 1) +
  geom_point(data = worldgroup_dat %>% filter(!species %in% "Natator_depressus"), aes(x = protected_perc, y = mean_mass, color = species), size = 3) +
  xlab("Nesting sites within protected area (%)") + ylab("Mean body mass (kg)") +
  scale_colour_manual(values = mycol) +
  facet_wrap(~ species, ncol = 2, scales = 'free') +
  ggrepel::geom_text_repel(data = worldgroup_dat %>% filter(!species %in% "Natator_depressus"), aes(x = protected_perc, y = mean_mass, label = name), size = 2, colour = "black") +
  mytheme + theme(legend.position = "bottom")

cowplot::plot_grid(spatial_plot, protect_plot, ncol = 1, rel_heights = c(0.6, 1))

