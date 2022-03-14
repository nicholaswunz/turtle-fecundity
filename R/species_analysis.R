# SET UP
# Install and load packages
library(ggplot2)
library(cowplot)
library(dplyr) 
library(tidyr)
library(reshape2)
library(brms)
library(rstan)
library(bayesplot)
library(rotl)
library(ape)
library(ggtree)

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(1, 1, 1, 1), units = , "cm"), 
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

mycol <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

# Set directory
setwd('repository space')

# Load dataset
meta_data <- read.csv("meta_raw_data.csv")
egg_data  <- read.csv("egg_energy_data.csv")
phylo_cor <- readRDS("phylo_cor.rds")
phyo_tree <- readRDS("phyo_tree.rds")

# Clean raw data
clean_meta_data <- meta_data %>% 
  dplyr::select(study_ID, species, red_list, population_level, pop_ID, country, lat, lon, temp,
                egg_volume_cm3, hatchling_length_mm, pred_clutch, pred_mass_kg) %>% 
  dplyr::mutate(lnMass   = log(pred_mass_kg),
                lnClutch = log(pred_clutch),
                lnVolume = log(egg_volume_cm3),
                lnHatch  = log(hatchling_length_mm),
                species  = as.character(species),
                animal   = as.factor(species)) # rename animal for phylogenetic random effects, vs within species effect

clean_egg_data  <- egg_data %>% 
  dplyr::select(study_ID, species, country, energy_kJ, egg_mass_g) %>% 
  dplyr::mutate(lnEnergy = log(energy_kJ),
                lnMass   = log(egg_mass_g))

species         <- sort(unique(clean_meta_data$species)) # generating list of species

# Filter by within and between populations
clean_meta_within  <- clean_meta_data %>% 
  dplyr::filter(population_level == "within population")

clean_meta_between <- clean_meta_data %>% 
  dplyr::filter(population_level == "between population")

## SUMMARY DATASET ## -------------------------------------------------------------------
nrow(clean_meta_data)
unique(clean_meta_data$study_ID)

meta_data %>% 
  dplyr::filter(population_level == "between population") %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(adult_length_cm     = length(adult_length_cm[!is.na(adult_length_cm)]),
                   adult_mass_kg       = length(adult_mass_kg[!is.na(adult_mass_kg)]),
                   egg_diameter_mm     = length(egg_diameter_mm[!is.na(egg_diameter_mm)]),
                   egg_mass_g          = length(egg_mass_g[!is.na(egg_mass_g)]),
                   hatchling_length_mm = length(hatchling_length_mm[!is.na(hatchling_length_mm)]),
                   hatchling_mass_g    = length(hatchling_mass_g[!is.na(hatchling_mass_g)]),
                   clutch_size         = length(clutch_size[!is.na(clutch_size)]),
                   clutch_mass_g       = length(clutch_mass_g[!is.na(clutch_mass_g)]),
                   clutch_freq         = length(clutch_freq[!is.na(clutch_freq)]))

country <- as.data.frame(clean_meta_between %>%
                           dplyr::group_by(species, country) %>% 
                           dplyr::summarise(mean_temp = mean(temp[!is.na(temp)]),
                                            mass_mean = mean(pred_mass_kg[!is.na(pred_mass_kg)]),
                                            mass_sd   = sd(pred_mass_kg[!is.na(pred_mass_kg)]),
                                            mass_n    = length(pred_mass_kg[!is.na(pred_mass_kg)])))

## ACROSS SPECIES ANALYSIS ## ------------------------------------------------------
# General STAN spec
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores available to use

# Set priors
meta_priors <- c(prior(normal(0, 10), "b"), # mean of 0 and SD of 10 (wide distribution)
                 prior(normal(0, 10), "Intercept"), # mean of 0 and SD of 10 (if log then wide negative disribution becomes infinite)
                 prior(student_t(3, 0, 20), "sd"), # class of random effect deviation to calcuate - has to be postiive (no negative SD)
                 prior(student_t(3, 0, 20), "sigma")) # residual SD parameter

# Fecundity-mass model
set.seed(10)
between_fecundity_model <- brms::brm(lnClutch ~ lnMass + (1 + lnMass | species) + (1 | gr(animal, cov = phylo)), 
                                    data    = clean_meta_between, 
                                    family  = gaussian(), 
                                    data2   = list(phylo = phylo_cor),
                                    prior   = meta_priors,
                                    chains  = 4, 
                                    cores   = 4, 
                                    iter    = 5e3, 
                                    warmup  = 2.5e3, 
                                    control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(between_fecundity_model)

# Estimate phylogenetic signal (delta) via hypothesis method
(hyp <- hypothesis(between_fecundity_model, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL))

# Check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
clean_meta_between2 <- clean_meta_between %>% tidyr::drop_na(lnClutch, lnMass)
clutch_y            <- clean_meta_between2$lnClutch
clutch_yrep         <- posterior_predict(between_fecundity_model, draws = 500)
species_pred        <- clean_meta_between2$species
bayesplot::ppc_scatter_avg(clutch_y, clutch_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(4.6, 6.7), y = c(4.6, 6.7)) + geom_point(aes(colour = species_pred), size = 2.5) + scale_colour_manual(values = mycol)

# Extract the full posterior distribution of species-level intercepts and slopes (1e4 posterior samples)
between_fecundity_post <- brms::posterior_samples(between_fecundity_model)
head(between_fecundity_post) # b_* indicates the fixed effects; r_animal* intercept-level phylogenetic random effects; r_species* intercept- and slope-level within-species random effects;

# for each row in the posterior matrix, obtain the "animal" random effect estimates for the intercept as well as the "species" random effects on the intercept and slope
# intercept = random_intercept_animal + random_intercept_species + fixed_intercept
# slope = random_slope_species + fixed_slope
out_list  <-  vector(mode = 'list', length = length(species))
for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(between_fecundity_post), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(between_fecundity_post)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(between_fecundity_post[i, 'b_Intercept'] + between_fecundity_post[i, paste0('r_animal[', species[j], ',Intercept]')] + between_fecundity_post[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  between_fecundity_post[i, 'b_lnMass'] + between_fecundity_post[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}

fecundity_out_data <- do.call('cbind.data.frame', out_list)
fecundity_post     <- sample_n(fecundity_out_data, 4000) # Extract 4000 random rows for reproductive output analysis

# EggVol-mass model
between_eggvol_model <- brms::brm(lnVolume ~ lnMass + (1 + lnMass | species) + (1 | gr(animal, cov = phylo)), 
                                 data    = clean_meta_between, 
                                 family  = gaussian(), 
                                 data2   = list(phylo = phylo_cor),
                                 prior   = meta_priors,
                                 chains  = 4, 
                                 cores   = 4, 
                                 iter    = 5e3, 
                                 warmup  = 2.5e3, 
                                 control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(between_eggvol_model)

# Estimate phylogenetic signal (delta) via hypothesis method
(hyp <- hypothesis(between_eggvol_model, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL))

# Check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
clean_meta_between2 <- clean_meta_between %>% tidyr::drop_na(lnVolume, lnMass)
egg_y               <- clean_meta_between2$lnVolume
egg_yrep            <- posterior_predict(between_eggvol_model, draws = 500)
species_pred        <- clean_meta_between2$species
bayesplot::ppc_scatter_avg(egg_y, egg_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(3, 4.7), y = c(3, 4.7)) + geom_point(aes(colour = species_pred), size = 2.5) + scale_colour_manual(values = mycol)

# Extract the full posterior distribution of species-level intercepts and slopes (1e4 posterior samples)
between_egg_post <- brms::posterior_samples(between_eggvol_model)
out_list         <-  vector(mode = 'list', length = length(species))
for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(between_egg_post), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(between_egg_post)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(between_egg_post[i, 'b_Intercept'] + between_egg_post[i, paste0('r_animal[', species[j], ',Intercept]')] + between_egg_post[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  between_egg_post[i, 'b_lnMass'] + between_egg_post[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}
egg_out_data <-  do.call('cbind.data.frame', out_list)

# Extract 4000 random rows for reproductive output analysis
egg_vol_post <- sample_n(egg_out_data, 4000)

# Combine fecundity_post and egg_post together
fecundity_int <- fecundity_post %>%
  reshape2::melt(measure.vars  = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas", "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata", "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea", "Intercept_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "fecundity_int") %>%
  dplyr::mutate(species = as.factor(substr(species, 11, 100))) %>% # remove the word "intercept"
  dplyr::select(species, fecundity_int)
  
fecundity_slope <- fecundity_post %>%
  reshape2::melt(measure.vars  = c("Slope_Caretta_caretta","Slope_Chelonia_mydas", "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata", "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea", "Slope_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "fecundity_slope") %>%
  dplyr::mutate(species = as.factor(substr(species, 7, 100))) %>% # remove the word "slope_"
  dplyr::select(fecundity_slope)

eggvolume_int <- egg_vol_post %>%
  reshape2::melt(measure.vars  = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas", "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata", "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea", "Intercept_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "egg_vol_int") %>%
  dplyr::mutate(species = as.factor(substr(species, 11, 100))) %>% # remove the word "intercept"
  dplyr::select(egg_vol_int)

eggvolume_slope <- egg_vol_post %>%
  reshape2::melt(measure.vars  = c("Slope_Caretta_caretta","Slope_Chelonia_mydas", "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata", "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea", "Slope_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "egg_vol_slope") %>%
  dplyr::mutate(species = as.factor(substr(species, 7, 100))) %>% # remove the word "slope_"
  dplyr::select(egg_vol_slope)

# Merge fecundity and egg volume data, calculate total reproductive output
slope_post <- data.frame(fecundity_int, fecundity_slope, eggvolume_int, eggvolume_slope) %>%
  dplyr::mutate(clutchvolume_int   = fecundity_int * egg_vol_int,
                clutchvolume_slope = fecundity_slope + egg_vol_slope)

# Egg-energy model
eggenergy_model <- brm(lnEnergy ~ lnMass + (1 | species), 
                       data    = clean_egg_data,  
                       family  = gaussian(),
                       prior   = meta_priors,
                       chains  = 4, 
                       cores   = 4, 
                       iter    = 5e3, 
                       warmup  = 2.5e3, 
                       control = list(adapt_delta = 0.999, max_treedepth = 20))

summary(eggenergy_model)

# Check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
energy_y            <- clean_egg_data$lnEnergy
energy_yrep         <- posterior_predict(eggenergy_model, draws = 500)
species_pred        <- clean_egg_data$species
bayesplot::ppc_scatter_avg(energy_y, energy_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(4.6, 6.7), y = c(4.6, 6.7)) + geom_point(aes(colour = species_pred), size = 2.5) + scale_colour_manual(values = mycol)

# Calculate total reproductive-energy output
slope_post <- slope_post %>%
  dplyr::mutate(egg_energy_int     = exp(brms::fixef(eggenergy_model)['Intercept', 'Estimate']),
                egg_energy_slope   = brms::fixef(eggenergy_model)['lnMass', 'Estimate'],
                clutchenergy_int   = fecundity_int * egg_energy_int * (egg_vol_int ^ egg_energy_slope),
                clutchenergy_slope = fecundity_slope + (egg_vol_slope * egg_energy_slope))

# Summarise output
energy_sum <- as.data.frame(slope_post %>%
  dplyr::group_by(species) %>% 
  dplyr::summarise(mean_fecundity_int      = mean(log(fecundity_int)),
                   mean_fecundity_slope    = mean(fecundity_slope),
                   mean_egg_vol_int        = mean(log(egg_vol_int)),
                   mean_egg_vol_slope      = mean(egg_vol_slope),
                   mean_clutchvol_int      = mean(clutchvolume_int),
                   mean_clutchvol_slope    = mean(clutchvolume_slope),
                   ci_clutchvol_slope      = bayestestR::ci(clutchvolume_slope),
                   mean_clutchenergy_int   = mean(clutchenergy_int),
                   mean_clutchenergy_slope = mean(clutchenergy_slope),
                   ci_clutchenergy_slope   = bayestestR::ci(clutchenergy_slope)))

# Create matrix of min and max values per group
mass_range <- clean_meta_between %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(min = min(lnMass[!is.na(lnMass)]),
                   max = max(lnMass[!is.na(lnMass)])) %>%
  as.data.frame()

fecundity_marg_eff <- as.data.frame(ggeffects::ggpredict(between_fecundity_model, terms = "lnMass"))
eggvol_marg_eff    <- as.data.frame(ggeffects::ggpredict(between_eggvol_model, terms = "lnMass"))

# Predict mass within min-max values
pred_df <- as.data.frame(mass_range %>%
                           mutate(mass_pred = purrr::map2(min, max, seq, length.out = 20)) %>%
                           unnest(cols = mass_pred) %>%
                           dplyr::mutate(mean_fecundity_int   = energy_sum$mean_fecundity_int[match(species, energy_sum$species)],
                                         mean_fecundity_slope = energy_sum$mean_fecundity_slope[match(species, energy_sum$species)],
                                         mean_egg_vol_int     = energy_sum$mean_egg_vol_int[match(species, energy_sum$species)],
                                         mean_egg_vol_slope   = energy_sum$mean_egg_vol_slope[match(species, energy_sum$species)]) %>%
                           dplyr::mutate(clutch_pred          = mean_fecundity_int + mean_fecundity_slope * mass_pred,
                                         eggvol_pred          = mean_egg_vol_int + mean_egg_vol_slope * mass_pred)) 
# Fig 1a
clutch_plot <- ggplot() +
  geom_ribbon(data = fecundity_marg_eff, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = fecundity_marg_eff, aes(x = x, y = conf.low), linetype = "dashed", alpha = 0.6) +
  geom_line(data = fecundity_marg_eff, aes(x = x, y = conf.high), linetype = "dashed", alpha = 0.6) +
  geom_line(data = fecundity_marg_eff, aes(x = x, y = predicted), size = 1, alpha = 0.6) +
  geom_point(data = clean_meta_between, aes(x = lnMass, y = lnClutch, fill = species, colour = species), alpha = 0.5, show.legend = FALSE) +
  geom_line(data = pred_df, aes(x = mass_pred, y = clutch_pred, colour = species), size = 1.5) +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab(expression("Fecundity (total clutch output)")) +
  scale_x_continuous(breaks = c(3.401197, 3.912023, 4.60517, 5.298317,  6.214608), # to get 30, 50, 100, 200, 500 kg
                     labels = round(c(exp(3.401197),exp(3.912023), exp(4.60517), exp(5.298317), exp(6.214608)), digits = 1)) +
  scale_y_continuous(breaks = c(0, 4.60517, 5.298317, 6.214608, 6.907755),# to get 1, 100, 200, 500, 1000 n of eggs
                     labels = round(c(exp(0), exp(4.60517), exp(5.298317), exp(6.214608), exp(6.907755)), digits = 1)) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol) +
  mytheme + theme(legend.position = "top")

# Fig. 1b
eggvol_plot <- ggplot() +
  geom_ribbon(data = eggvol_marg_eff, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = eggvol_marg_eff, aes(x = x, y = conf.low), linetype = "dashed", alpha = 0.6) +
  geom_line(data = eggvol_marg_eff, aes(x = x, y = conf.high), linetype = "dashed", alpha = 0.6) +
  geom_line(data = eggvol_marg_eff, aes(x = x, y = predicted), size = 1, alpha = 0.6) +
  geom_point(data = clean_meta_between, aes(x = lnMass, y = lnVolume, fill = species, colour = species), alpha = 0.5, show.legend = FALSE) +
  geom_line(data = pred_df, aes(x = mass_pred, y = eggvol_pred, colour = species), size = 1.5) +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab(expression("Egg volume (cm"^"3"*")")) +
  scale_x_continuous(breaks = c(3.401197, 3.912023, 4.60517, 5.298317,  6.214608), # to get 30, 50, 100, 200, 500 kg
                     labels = round(c(exp(3.401197),exp(3.912023), exp(4.60517), exp(5.298317), exp(6.214608)), digits = 1)) +
  scale_y_continuous(breaks = c(2.302585, 3.218876, 3.912023, 4.317488, 4.60517),# to get 10, 25, 50, 75, 100 cm3
                     labels = round(c(exp(2.302585), exp(3.218876), exp(3.912023), exp(4.317488), exp(4.60517)), 1)) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol) +
  mytheme + theme(legend.position = "top")

between_mass_plot <- plot_grid(clutch_plot, eggvol_plot, align = "h", axis = "bt")

# Plot Fig 1c and d
energy_sum2    <- tibble::remove_rownames(energy_sum) %>% tibble::column_to_rownames(var = "species")
name_match     <- intersect(phyo_tree$tip.label, rownames(energy_sum2)) # matching tree tip label with rownames
name_corrected <- energy_sum2[name_match,] # attach intersect with dataframe
tree_plot      <- ggtree::ggtree(phyo_tree, size = 1) + geom_tiplab(size = 2, hjust = -0.1)

clutchvol_slope    <- data.frame(species = phyo_tree$tip.label, 
                                 mean = name_corrected$mean_clutchvol_slope, 
                                 xmin = name_corrected$ci_clutchvol_slope$CI_low, 
                                 xmax = name_corrected$ci_clutchvol_slope$CI_high)

clutchenergy_slope <- data.frame(species = phyo_tree$tip.label, 
                                 mean = name_corrected$mean_clutchenergy_slope, 
                                 xmin = name_corrected$ci_clutchenergy_slope$CI_low, 
                                 xmax = name_corrected$ci_clutchenergy_slope$CI_high)

clutchvol_plot     <- facet_plot(tree_plot + xlim_tree(1.3), 
                                 panel = "Reproductive-output exponents",
                                 data = clutchvol_slope, 
                                 geom = ggstance::geom_pointrangeh, 
                                 mapping = aes(x = mean, xmin = xmin, xmax = xmax, colour = species), size = 0.7) +
  theme_tree2() + geom_vline(xintercept = 1, linetype = "dashed") + scale_colour_manual(values = mycol) + mytheme

phylo_plot <- facet_plot(clutchvol_plot + xlim_tree(1.3), 
                         panel = "Reproductive-energy exponents",
                         data = clutchenergy_slope, 
                         geom = ggstance::geom_pointrangeh, 
                         mapping = aes(x = mean, xmin = xmin, xmax = xmax, colour = species), size = 0.7)

# Fig 1
cowplot::plot_grid(between_mass_plot, phylo_plot, ncol = 1, align = "h", axis = "bt", rel_heights = c(1, 0.6))


## HATCHLING ANALYSIS ## -------------------------------------------------------------------------
# Hatchling-mass model
hatchling_model <- brms::brm(lnHatch ~ lnMass + (1 + lnMass | species) + (1 | gr(animal, cov = phylo)), 
                                  data    = clean_meta_between, 
                                  family  = gaussian(), 
                                  data2   = list(phylo = phylo_cor),
                                  prior   = meta_priors,
                                  chains  = 4, 
                                  cores   = 4, 
                                  iter    = 5e3, 
                                  warmup  = 2.5e3, 
                                  control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(hatchling_model)

# Estimate phylogenetic signal (delta) via hypothesis method
(hyp <- hypothesis(hatchling_model, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL))

# Check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
clean_meta_between2 <- clean_meta_between %>% tidyr::drop_na(lnHatch, lnMass)
hatch_y             <- clean_meta_between2$lnHatch
hatch_yrep          <- posterior_predict(hatchling_model, draws = 500)
species_pred        <- clean_meta_between2$species
bayesplot::ppc_scatter_avg(hatch_y, hatch_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(3.6, 4.2), y = c(3.6, 4.2)) + geom_point(aes(colour = species_pred), size = 2.5) + scale_colour_manual(values = mycol)

# Extract the full posterior distribution of species-level intercepts and slopes (1e4 posterior samples)
hatchling_post <- brms::posterior_samples(hatchling_model)
out_list       <-  vector(mode = 'list', length = length(species))
for (j in seq_along(species)) {
  out_list[[j]]  <-  data.frame(matrix(NA, nrow(hatchling_post), 2))
  dat_names      <-  paste0(c('Intercept_', 'Slope_'), species[j])
  names(out_list[[j]])  <-  dat_names
  for (i in 1:nrow(hatchling_post)) {
    out_list[[j]][i, dat_names[1]]  <-  exp(hatchling_post[i, 'b_Intercept'] + hatchling_post[i, paste0('r_animal[', species[j], ',Intercept]')] + hatchling_post[i, paste0('r_species[', species[j], ',Intercept]')])
    out_list[[j]][i, dat_names[2]]  <-  hatchling_post[i, 'b_lnMass'] + hatchling_post[i, paste0('r_species[', species[j], ',lnMass]')]
  }
}
hatchling_out_data <-  do.call('cbind.data.frame', out_list)

# Extract 4000 random rows for reproductive output analysis
hatchling_post <- sample_n(hatchling_out_data, 4000)

# Combine fecundity_post and egg_post together
hatchling_int <- hatchling_post %>%
  reshape2::melt(measure.vars  = c("Intercept_Caretta_caretta","Intercept_Chelonia_mydas", "Intercept_Dermochelys_coriacea", "Intercept_Eretmochelys_imbricata", "Intercept_Lepidochelys_kempii", "Intercept_Lepidochelys_olivacea", "Intercept_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "hatchling_int") %>%
  dplyr::mutate(species = as.factor(substr(species, 11, 100))) %>% # remove the word "intercept"
  dplyr::select(species, hatchling_int)

hatchling_slope <- hatchling_post %>%
  reshape2::melt(measure.vars  = c("Slope_Caretta_caretta","Slope_Chelonia_mydas", "Slope_Dermochelys_coriacea", "Slope_Eretmochelys_imbricata", "Slope_Lepidochelys_kempii", "Slope_Lepidochelys_olivacea", "Slope_Natator_depressus"),
                 variable.name = "species", 
                 value.name    = "hatchling_slope") %>%
  dplyr::mutate(species = as.factor(substr(species, 7, 100))) %>% # remove the word "slope_"
  dplyr::select(hatchling_slope)

# Merge fecundity and egg volume data, calculate total reproductive output
hatchling_group_post <- data.frame(hatchling_int, hatchling_slope)

# Summarise output
hatchling_sum <- as.data.frame(hatchling_group_post %>%
                              dplyr::group_by(species) %>% 
                              dplyr::summarise(mean_hatchling_int   = mean(log(hatchling_int)),
                                               mean_hatchling_slope = mean(hatchling_slope)))

# Create matrix of min and max values per group
mass_range <- clean_meta_between %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(min = min(lnMass[!is.na(lnMass)]),
                   max = max(lnMass[!is.na(lnMass)])) %>%
  as.data.frame()

hatchling_marg_eff <- as.data.frame(ggeffects::ggpredict(hatchling_model, terms = "lnMass"))

# Predict mass within min-max values
pred_df <- as.data.frame(mass_range %>%
                           mutate(mass_pred = purrr::map2(min, max, seq, length.out = 20)) %>%
                           unnest(cols = mass_pred) %>%
                           dplyr::mutate(mean_hatchling_int   = hatchling_sum$mean_hatchling_int[match(species, hatchling_sum$species)],
                                         mean_hatchling_slope = hatchling_sum$mean_hatchling_slope[match(species, hatchling_sum$species)]) %>%
                           dplyr::mutate(hatchling_pred          = mean_hatchling_int + mean_hatchling_slope * mass_pred))

# Fig S6
ggplot() +
  geom_ribbon(data = hatchling_marg_eff, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.1) +
  geom_line(data = hatchling_marg_eff, aes(x = x, y = conf.low), linetype = "dashed", alpha = 0.6) +
  geom_line(data = hatchling_marg_eff, aes(x = x, y = conf.high), linetype = "dashed", alpha = 0.6) +
  geom_line(data = hatchling_marg_eff, aes(x = x, y = predicted), size = 1, alpha = 0.6) +
  geom_point(data = clean_meta_between, aes(x = lnMass, y = lnHatch, fill = species, colour = species), alpha = 0.5, show.legend = FALSE) +
  geom_line(data = pred_df, aes(x = mass_pred, y = hatchling_pred, colour = species), size = 1.5) +
  xlab(expression("Female body mass (kg)")) + 
  ylab(expression("Hatchling length (mm)")) +
  scale_x_continuous(breaks = c(3.401197, 3.912023, 4.60517, 5.298317,  6.214608), # to get 30, 50, 100, 200, 500 kg
                     labels = round(c(exp(3.401197),exp(3.912023), exp(4.60517), exp(5.298317), exp(6.214608)), digits = 1)) +
  scale_y_continuous(breaks = c(3.401197, 3.688879, 3.912023, 4.094345, 4.248495),# to get 30, 40, 50, 60, 70 mm
                     labels = round(c(exp(3.401197), exp(3.688879), exp(3.912023), exp(4.094345), exp(4.248495)), digits = 1)) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol) +
  mytheme + theme(legend.position = "top")


## WITHIN-POPULATION ANALYSIS ## -------------------------------------------------------------------------
# Data summary
as.data.frame(clean_meta_within %>%
  dplyr::group_by(species, country) %>%
  dplyr::summarise(clutch_n = length(pred_clutch[!is.na(pred_clutch)]),
                   eggvol_n = length(egg_volume_cm3[!is.na(egg_volume_cm3)])))

# Fecundity-mass model 
within_fecundity_model <- brms::brm(lnClutch ~ lnMass + (1 + lnMass | species / pop_ID) + (1 | gr(animal, cov = phylo)), 
                                    data    = clean_meta_within, 
                                    family  = gaussian(), 
                                    data2   = list(phylo = phylo_cor),
                                    prior   = meta_priors,
                                    chains  = 4, 
                                    cores   = 4, 
                                    iter    = 5e3, 
                                    warmup  = 2.5e3, 
                                    control = list(adapt_delta = 0.999, max_treedepth = 15))
summary(within_fecundity_model)

# Estimate phylogenetic signal (delta) via hypothesis method
(hyp <- hypothesis(within_fecundity_model, "sd_animal__Intercept^2 / (sd_animal__Intercept^2 + sigma^2) = 0", class = NULL))

fixefs_within_fecundity <- fixef(within_fecundity_model) # fixed effects summary (all data)
clutch_animal           <- as.data.frame(ranef(within_fecundity_model)$animal[, 'Estimate', 'Intercept']) %>% dplyr::rename(animal_intercept = 'ranef(within_fecundity_model)$animal[, "Estimate", "Intercept"]')
clutch_pop              <- as.data.frame(ranef(within_fecundity_model)$"species:pop_ID"[, 'Estimate', 'Intercept']) %>% dplyr::rename(estimate_intercept = 'ranef(within_fecundity_model)$"species:pop_ID"[, "Estimate", "Intercept"]')

clutch_intercept <- as.data.frame(ranef(within_fecundity_model)$"species:pop_ID"[, 'Estimate', 'Intercept']) %>%
  dplyr::rename(estimate_intercept = 'ranef(within_fecundity_model)$"species:pop_ID"[, "Estimate", "Intercept"]') %>%
  dplyr::mutate(ID = row.names(clutch_pop)) %>%
  tidyr::separate(ID, into = c("species", "population"), sep = "_(?=[^_]*_[^_]+$)") %>% 
  dplyr::mutate(animal_estimate = clutch_animal$animal_intercept[match(species, row.names(clutch_animal))],
                fecundity_int   = exp(fixefs_within_fecundity['Intercept', 'Estimate'] + animal_estimate + estimate_intercept)) #intercept = random_intercept_animal + random_intercept_species + fixed_intercept

clutch_slope <- as.data.frame(ranef(within_fecundity_model)$"species:pop_ID"[, , 'lnMass']) %>%
  dplyr::rename(estimate_lnMass = Estimate,
                Q2.5_lnMass     = Q2.5,
                Q97.5_lnMass    = Q97.5) %>%
  dplyr::mutate(ID = row.names(clutch_pop)) %>%
  dplyr::mutate(fecundity_slope_est = fixefs_within_fecundity['lnMass', 'Estimate'] + estimate_lnMass,
                fecundity_slope_low = fixefs_within_fecundity['lnMass', 'Estimate'] + Q2.5_lnMass,
                fecundity_slope_upp = fixefs_within_fecundity['lnMass', 'Estimate'] + Q97.5_lnMass) # Slope = random_slope_species + fixed_slope

# EggVol-mass model (only 5 CC pop, 3 CM pop, and 1 DC pop available) - not appropriate with phylogenetically-corrected model
# Non-phylo analysis
within_eggvol_model <- brms::brm(lnVolume ~ lnMass + (1 + lnMass | species / pop_ID), 
                                 data    = clean_meta_within, 
                                 family  = gaussian(), 
                                 prior   = meta_priors,
                                 chains  = 4, 
                                 cores   = 4, 
                                 iter    = 5e3, 
                                 warmup  = 2.5e3, 
                                 control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(within_eggvol_model)

fixefs_within_eggvol <- fixef(within_eggvol_model) # fixed effects summary (all data)
effvol_pop           <- as.data.frame(ranef(within_eggvol_model)$"species:pop_ID"[, 'Estimate', 'Intercept']) %>% dplyr::rename(estimate_intercept = 'ranef(within_eggvol_model)$"species:pop_ID"[, "Estimate", "Intercept"]')

eggvol_intercept <- as.data.frame(ranef(within_eggvol_model)$"species:pop_ID"[, 'Estimate', 'Intercept']) %>%
  dplyr::rename(estimate_intercept = 'ranef(within_eggvol_model)$"species:pop_ID"[, "Estimate", "Intercept"]') %>%
  dplyr::mutate(ID = row.names(effvol_pop)) %>%
  tidyr::separate(ID, into = c("species", "population"), sep = "_(?=[^_]*_[^_]+$)") %>% 
  dplyr::mutate(egg_vol_int = exp(fixefs_within_eggvol['Intercept', 'Estimate'] + estimate_intercept)) #intercept = random_intercept_species + fixed_intercept

eggvol_slope <- as.data.frame(ranef(within_eggvol_model)$"species:pop_ID"[, , 'lnMass']) %>%
  dplyr::rename(estimate_lnMass = Estimate,
                Q2.5_lnMass     = Q2.5,
                Q97.5_lnMass    = Q97.5) %>%
  dplyr::mutate(ID = row.names(effvol_pop)) %>%
  dplyr::mutate(egg_vol_slope_est = fixefs_within_eggvol['lnMass', 'Estimate'] + estimate_lnMass,
                egg_vol_slope_low = fixefs_within_eggvol['lnMass', 'Estimate'] + Q2.5_lnMass,
                egg_vol_slope_upp = fixefs_within_eggvol['lnMass', 'Estimate'] + Q97.5_lnMass)  # Slope = random_slope_species + fixed_slope

pop_grouped <- merge(clutch_intercept, clutch_slope, by = "row.names") %>%
  dplyr::mutate(egg_vol_int       = eggvol_intercept$egg_vol_int[match(ID, row.names(eggvol_intercept))],
                egg_vol_slope_est = eggvol_slope$egg_vol_slope_est[match(ID, row.names(eggvol_slope))],
                egg_vol_slope_low = eggvol_slope$egg_vol_slope_low[match(ID, row.names(eggvol_slope))],
                egg_vol_slope_upp = eggvol_slope$egg_vol_slope_upp[match(ID, row.names(eggvol_slope))]) %>%
  dplyr::select(ID, species, population, fecundity_int, fecundity_slope_est, fecundity_slope_low, fecundity_slope_upp, 
                egg_vol_int, egg_vol_slope_est, egg_vol_slope_low, egg_vol_slope_upp) %>% 
  dplyr::mutate(clutchvolume_int   = fecundity_int * egg_vol_int,
                clutchvolume_slope_est = fecundity_slope_est + egg_vol_slope_est,
                clutchvolume_slope_low = fecundity_slope_est + egg_vol_slope_low,
                clutchvolume_slope_upp = fecundity_slope_est + egg_vol_slope_upp)

ggplot(eggvol_slope, aes(x = ID, y = egg_vol_slope_est)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = egg_vol_slope_low, ymax = egg_vol_slope_upp)) +
  xlab(NULL) +
  coord_flip() +
  mytheme()

## PLOT POSTERIOR DISTRIBUTION BY SPECIES AND POPULATION ## ----------------------------------------------------------------------------------------
# Fecundity
within_fecundity_post <- posterior_samples(within_fecundity_model) # extract posterior
clean_meta_within     <- clean_meta_within %>% tidyr::unite("pop_unique", c("species", "pop_ID"), sep = "_", remove = FALSE) # create pop_unique column
clean_fecund_within   <- subset(clean_meta_within, clean_meta_within$lnClutch != "NA") # remove fecundity data with NA

# Loop to extract within population-level posterior of slope
pop_unique            <- sort(unique(clean_fecund_within$pop_unique)) # pop ID names
out_list              <- vector(mode = 'list', length = length(pop_unique))
for (j in seq_along(pop_unique)) {
  out_list[[j]]                 <- data.frame(matrix(NA, nrow(within_fecundity_post), 1))
  dat_names                     <-  paste0(pop_unique[j])
  names(out_list[[j]])          <-  dat_names
  out_list[[j]][, dat_names[1]] <-  within_fecundity_post[, paste0('r_species:pop_ID[', pop_unique[j], ',lnMass]')]
}

pop_fecundity_slope_out  <- do.call('cbind.data.frame', out_list)
pop_fecundity_slope_post <- sample_n(pop_fecundity_slope_out, 4000) # Extract 4000 random rows for reproductive output analysis

# Loop to extract within species-level posterior of slope
species_unique <- sort(unique(clean_fecund_within$species)) # species ID names
out_list       <-  vector(mode = 'list', length = length(species_unique))
for (j in seq_along(species_unique)) {
  out_list[[j]]        <- data.frame(matrix(NA, nrow(within_fecundity_post), 1))
  dat_names            <- paste0(species_unique[j])
  names(out_list[[j]]) <- dat_names
  out_list[[j]][, dat_names[1]] <- within_fecundity_post[, 'b_lnMass'] + within_fecundity_post[, paste0('r_species[', species_unique[j], ',lnMass]')]
}

sp_fecundity_slope_out  <- do.call('cbind.data.frame', out_list)
sp_fecundity_slope_post <- sample_n(sp_fecundity_slope_out, 4000) # Extract 4000 random rows for reproductive output analysis

# Tidy data
pop_fecundity_slope_post_long <- pop_fecundity_slope_post %>%
  tidyr::pivot_longer(cols = everything(),names_to = "pop_unique", values_to = "slope") %>%
  dplyr::filter(pop_unique != "Caretta_caretta_")
pop_fecundity_slope_post_long$pop_ID  <- clean_fecund_within$pop_ID[match(pop_fecundity_slope_post_long$pop_unique, clean_fecund_within$pop_unique)]
pop_fecundity_slope_post_long$species <- clean_fecund_within$species[match(pop_fecundity_slope_post_long$pop_unique, clean_fecund_within$pop_unique)]

sp_fecundity_slope_post_long <- sp_fecundity_slope_post %>%
  tidyr::pivot_longer(cols = everything(),names_to = "species", values_to = "slope")

# Plot Fig 7a fecundity slopes
pop_fecundity_slope_post_long %>%
  ggplot() +
  geom_density(data = sp_fecundity_slope_post_long, aes(slope), fill = "black", colour = NA, alpha = 0.3) +
  geom_density(aes(slope, group = pop_ID, fill = pop_ID, colour = pop_ID), alpha = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") + # isometric relationship
  xlim(0, 3) +
  xlab(expression("Size-fecundity scaling exponent"~(italic(beta)[1]))) +
  facet_grid(species ~ .) +
  mytheme() + theme(legend.position = "bottom")

# Egg volume
within_eggvol_post  <- posterior_samples(within_eggvol_model) # extract posterior
clean_eggvol_within <- subset(clean_meta_within, clean_meta_within$lnVolume != "NA") # remove egb volume data with NA

# Loop to extract within population level posterior of slope
pop_unique          <- sort(unique(clean_eggvol_within$pop_unique)) # pop ID names
out_list            <- vector(mode = 'list', length = length(pop_unique))
for (j in seq_along(pop_unique)) {
  out_list[[j]]                 <- data.frame(matrix(NA, nrow(within_eggvol_post), 1))
  dat_names                     <- paste0(pop_unique[j])
  names(out_list[[j]])          <- dat_names
  out_list[[j]][, dat_names[1]] <- within_eggvol_post[, 'b_lnMass'] + within_eggvol_post[, paste0('r_species:pop_ID[', pop_unique[j], ',lnMass]')]
}

pop_eggvol_slope_out  <- do.call('cbind.data.frame', out_list)
pop_eggvol_slope_post <- sample_n(pop_eggvol_slope_out, 4000) # Extract 4000 random rows for reproductive output analysis

# Loop to extract within species level posterior of slope
species_unique <- sort(unique(clean_eggvol_within$species)) # species ID names
out_list       <- vector(mode = 'list', length = length(species_unique))
for (j in seq_along(species_unique)) {
  out_list[[j]]                 <- data.frame(matrix(NA, nrow(within_fecundity_post), 1))
  dat_names                     <- paste0(species_unique[j])
  names(out_list[[j]])          <- dat_names
  out_list[[j]][, dat_names[1]] <- within_eggvol_post[, 'b_lnMass'] + within_eggvol_post[, paste0('r_species[', species_unique[j], ',lnMass]')]
}

sp_eggvol_slope_out  <- do.call('cbind.data.frame', out_list)
sp_eggvol_slope_post <- sample_n(sp_eggvol_slope_out, 4000) # Extract 4000 random rows for reproductive output analysis

# Tidy data
pop_eggvol_slope_post_long <- pop_eggvol_slope_post %>%
  tidyr::pivot_longer(cols = everything(),names_to = "pop_unique", values_to = "slope") %>%
  dplyr::filter(pop_unique != "Caretta_caretta_")
pop_eggvol_slope_post_long$pop_ID  <- clean_eggvol_within$pop_ID[match(pop_eggvol_slope_post_long$pop_unique, clean_eggvol_within$pop_unique)]
pop_eggvol_slope_post_long$species <- clean_eggvol_within$species[match(pop_eggvol_slope_post_long$pop_unique, clean_eggvol_within$pop_unique)]

sp_eggvol_slope_post_long <- sp_eggvol_slope_post %>%
  tidyr::pivot_longer(cols = everything(),names_to = "species", values_to = "slope")

# Plot Fig 7b egg volume slope
pop_eggvol_slope_post_long %>%
  ggplot() +
  geom_density(data = sp_eggvol_slope_post_long, aes(slope), fill = "black", colour = NA, alpha = 0.3) +
  geom_density(aes(slope, group = pop_ID, fill = pop_ID, colour = pop_ID), alpha = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") + # isometric relationship
  geom_vline(xintercept = 0) + # no relationship
  xlim(-1, 1) +
  xlab(expression("Size-egg volume scaling exponent"~(italic(beta)[1]))) +
  facet_grid(species ~ .) +
  mytheme() + theme(legend.position = "bottom")
