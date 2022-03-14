# SET UP
# Install and load packages from species_analysis.R

# Functions 
ss_tot <- function(x) {
  sum((x - mean(x))^2)
}

ss_res <- function(x, y) {
  sum((x - y)^2)
}

r_sqrd <- function(x, y) {
  1 - (ss_res(x, y) / ss_tot(x))
}

# Load data
nest_dat <- read.csv("nesting year.csv")

nest_dat <- nest_dat %>%
  dplyr::mutate(lnClutch = log(sum_eggs),
                lnMass   = log(mass_kg))

## DATA SUMMARY ## ------------------------------------------------------------------------------------
length(unique(nest_dat$ID)) # number of individuals tagged
nest.sum <- as.data.frame(nest_dat %>% 
  dplyr::group_by(year, recruit) %>% 
  dplyr::summarise(count  = length(unique(ID[!is.na(ID)])),
                   mass   = mean(mass_kg[!is.na(mass_kg)]),
                   clutch = round(mean(sum_eggs[!is.na(sum_eggs)]), digits = 0)))

# Calculate nesting error (how many one counts per year), then calculate the average frequency
nest_dat %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(one_nest = sum(count_year == 1),
                   total    = sum(count_year)) %>%
  dplyr::mutate(freq = one_nest / total * 100) %>%
  dplyr::summarise(mean = mean(freq), 
                   SD   = sd(freq))

## ESTIMATE TOTAL REPRODUCTIVE-ENERGY OUTPUT ## ------------------------------------------------------------------
# General STAN spec
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores available to use

# Set priors
quant_priors <- c(prior(normal(0, 10), "b"),
                  prior(normal(0, 10), "Intercept"),
                  prior(student_t(3, 0, 20), "sigma"))

# Fecundity models by quantile 
set.seed(10)
m_quantile_10     <- brm(bf(lnClutch ~ lnMass, quantile = 0.1), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))
m_quantile_50     <- brm(bf(lnClutch ~ lnMass, quantile = 0.5), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))
m_quantile_90     <- brm(bf(lnClutch ~ lnMass, quantile = 0.9), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))
m_quantile_10_iso <- brm(bf(sum_eggs ~ mass_kg, quantile = 0.1), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))
m_quantile_50_iso <- brm(bf(sum_eggs ~ mass_kg, quantile = 0.5), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))
m_quantile_90_iso <- brm(bf(sum_eggs ~ mass_kg, quantile = 0.9), data = nest_dat, family = asym_laplace(), prior = quant_priors, chains = 4, cores = 4, iter = 5e3, warmup = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 15))

bayes_R2(m_quantile_50)
bayes_R2(m_quantile_50_iso)

# Check observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
nest_dat2 <- nest_dat %>% drop_na(lnClutch)
clutch_y  <- nest_dat2$lnClutch
clutch_y2 <- nest_dat2$sum_eggs
yrep10    <- posterior_predict(m_quantile_10, draws = 500)
yrep50    <- posterior_predict(m_quantile_50, draws = 500)
yrep90    <- posterior_predict(m_quantile_90, draws = 500)
yrep50iso <- posterior_predict(m_quantile_50_iso, draws = 500)

pred10    <- bayesplot::ppc_scatter_avg(clutch_y, yrep10, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(3, 7.1), y = c(3, 7.1)) 
pred50    <- bayesplot::ppc_scatter_avg(clutch_y, yrep50, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(3, 7.1), y = c(3, 7.1)) 
pred90    <- bayesplot::ppc_scatter_avg(clutch_y, yrep90, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(3, 7.1), y = c(3, 7.1))
pred50iso <- bayesplot::ppc_scatter_avg(clutch_y2, yrep50iso, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(0, 1500), y = c(0, 1500)) 

cowplot::plot_grid(pred10, pred50, pred90, pred50iso, ncol = 2, align = "hv", axis = "tblr", labels = c("a", "b", "c", "d"))

# Calculate total reproductive-energy output
# Egg-energy-Volume model (KJ)
greenenergy_model <- brm(lnEnergy ~ lnMass, 
                       data    = clean_egg_data %>%
                         dplyr::filter(species == "Chelonia_mydas"),  
                       family  = gaussian(),
                       prior   = meta_priors,
                       chains  = 4, 
                       cores   = 4, 
                       iter    = 5e3, 
                       warmup  = 2.5e3, 
                       control = list(adapt_delta = 0.999, max_treedepth = 20))

summary(greenenergy_model)

exp_q10      <- exp(brms::fixef(m_quantile_10)["Intercept","Estimate"])
exp_q50      <- exp(brms::fixef(m_quantile_50)["Intercept","Estimate"])
exp_q90      <- exp(brms::fixef(m_quantile_90)["Intercept","Estimate"])
exp_q10_iso  <- brms::fixef(m_quantile_10_iso)["Intercept","Estimate"]
exp_q50_iso  <- brms::fixef(m_quantile_50_iso)["Intercept","Estimate"]
exp_q90_iso  <- brms::fixef(m_quantile_90_iso)["Intercept","Estimate"]
energy_post  <- brms::posterior_samples(greenenergy_model) # egg energy
quant10_post <- sample_n(brms::posterior_samples(m_quantile_10), 4000)
quant50_post <- sample_n(brms::posterior_samples(m_quantile_50), 4000)
quant90_post <- sample_n(brms::posterior_samples(m_quantile_90), 4000)
iso10_post   <- sample_n(brms::posterior_samples(m_quantile_10_iso), 4000)
iso50_post   <- sample_n(brms::posterior_samples(m_quantile_50_iso), 4000)
iso90_post   <- sample_n(brms::posterior_samples(m_quantile_90_iso), 4000)
energy_post  <- sample_n(energy_post, 4000)

str(slope_post) # from species_analysis.R

c_mydas_slope <- slope_post %>%
  dplyr::filter(species == "Chelonia_mydas") %>%
  dplyr::select(egg_vol_int, egg_vol_slope) %>%
  dplyr::mutate(energy_int           = exp(energy_post$b_Intercept),
                energy_slope         = energy_post$b_lnMass,
                fecundity10iso_int   = iso10_post$b_Intercept,
                fecundity10iso_slope = iso10_post$b_mass_kg,
                fecundity50iso_int   = iso50_post$b_Intercept,
                fecundity50iso_slope = iso50_post$b_mass_kg,
                fecundity90iso_int   = iso90_post$b_Intercept,
                fecundity90iso_slope = iso90_post$b_mass_kg,
                fecundity10_int      = exp(quant10_post$b_Intercept),
                fecundity10_slope    = quant10_post$b_lnMass,
                fecundity50_int      = exp(quant50_post$b_Intercept),
                fecundity50_slope    = quant50_post$b_lnMass,
                fecundity90_int      = exp(quant90_post$b_Intercept),
                fecundity90_slope    = quant90_post$b_lnMass,
                # calculate egg energy
                eggenergy_int       = energy_int * (egg_vol_int ^ energy_slope),
                eggenergy_slope     = egg_vol_slope * energy_slope)

# Calculate predicted isometry and allometry relationship for C. mydas (kJ to MJ)
mass_pred <- as.data.frame(seq(min(nest_dat$mass_kg - 10), max(nest_dat$mass_kg + 10))) %>% # in kg
  dplyr::rename(mass_pred = "seq(min(nest_dat$mass_kg - 10), max(nest_dat$mass_kg + 10))") %>%
  dplyr::mutate(hyper_10 = (mean(c_mydas_slope$fecundity10_int * c_mydas_slope$eggenergy_int) * mass_pred ^ mean(c_mydas_slope$fecundity10_slope + c_mydas_slope$eggenergy_slope) * 0.001),
                hyper_50 = (mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * mass_pred ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001),
                hyper_90 = (mean(c_mydas_slope$fecundity90_int * c_mydas_slope$eggenergy_int) * mass_pred ^ mean(c_mydas_slope$fecundity90_slope + c_mydas_slope$eggenergy_slope) * 0.001))

# Estimate reproductive-energy output based on known mass (For fig 3a)
(mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * 75 ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001) * 3 # for a small female (multiple by 3 to equal an 140 kg female)
(mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * 140 ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001) # for a large female 

# Fig 3a
green_hyper_plot <- ggplot(mass_pred, aes(x = mass_pred)) + mytheme() +
  geom_line(aes(y = hyper_10), size = 0.5, colour = "#348781", linetype = "dashed") +
  geom_line(aes(y = hyper_50), size = 1, colour = "#348781") +
  geom_line(aes(y = hyper_90), size = 0.5, colour = "#348781", linetype = "dashed") +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab(expression("Reproductive-energy output, "*italic(R)[i]*", (MJ)"))

# Compare prediction of iso and hyper with obs data
iso_10_post <- brms::posterior_epred(m_quantile_10_iso)
hyp_10_post <- brms::posterior_epred(m_quantile_10)
iso_50_post <- brms::posterior_epred(m_quantile_50_iso)
hyp_50_post <- brms::posterior_epred(m_quantile_50)
iso_90_post <- brms::posterior_epred(m_quantile_90_iso)
hyp_90_post <- brms::posterior_epred(m_quantile_90)

# Summarise fecundity_obs per year
obs_data <- nest_dat %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(fecundity_obs = sum(sum_eggs))

# Create empty matrix for loop
r2_iso_10 <- numeric(length = nrow(iso_10_post))
r2_hyp_10 <- numeric(length = nrow(hyp_10_post))
r2_iso_50 <- numeric(length = nrow(iso_50_post))
r2_hyp_50 <- numeric(length = nrow(hyp_50_post))
r2_iso_90 <- numeric(length = nrow(iso_90_post))
r2_hyp_90 <- numeric(length = nrow(hyp_90_post))

# create dataframe summarising predicted fecundity (iso and hyp) per year + r_sqrf with observed data
for (i in seq_len(nrow(iso_50_post))) {
  tmp <- data.frame(year = nest_dat$year,
                    iso_10  = iso_10_post[i, ],
                    hyp_10  = exp(hyp_10_post[i, ]),
                    iso_50  = iso_50_post[i, ],
                    hyp_50  = exp(hyp_50_post[i, ]),
                    iso_90  = iso_50_post[i, ],
                    hyp_90  = exp(hyp_50_post[i, ])) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(pred_iso_10 = sum(iso_10),
                     pred_hyp_10 = sum(hyp_10),
                     pred_iso_50 = sum(iso_50),
                     pred_hyp_50 = sum(hyp_50),
                     pred_iso_90 = sum(iso_90),
                     pred_hyp_90 = sum(hyp_90))
  r2_iso_50[i] <- r_sqrd(obs_data$fecundity_obs, tmp$pred_iso_50)
  r2_hyp_50[i] <- r_sqrd(obs_data$fecundity_obs, tmp$pred_hyp_50)
}

# Merge observed with predicted
pred_merged <- merge(obs_data, tmp, by = "year")

# Percent difference between hyper and iso models
min(pred_merged %>% dplyr::summarise(diff = (pred_hyp_50 - pred_iso_50) / pred_iso_50 * 100))
max(pred_merged %>% dplyr::summarise(diff = (pred_hyp_50 - pred_iso_50) / pred_iso_50 * 100))

# Fig 3b
obs_pred_plot <- pred_merged %>%
  ggplot(aes(x = fecundity_obs)) + 
  geom_abline(slope = 1, size = 1, linetype = "dashed") +
  geom_point(aes(y = pred_iso_50), size = 3, colour = "#858585", alpha = 0.5, shape = 17) +
  geom_point(aes(y = pred_hyp_50), size = 3, colour = "#348781", alpha = 0.5) +
  geom_smooth(aes(y = pred_iso_50), method = lm, se = FALSE, size = 1, colour = "#858585") +
  geom_smooth(aes(y = pred_hyp_50), method = lm, se = FALSE, size = 1, colour = "#348781") +
  lims(x = c(5000, 125000), y = c(5000, 125000)) +
  ylab("Predicted population fecundity") + xlab("Observed population fecundity") +
  mytheme

# R2
ggplot() +
  geom_density(data = as.data.frame(r2_iso), aes(x = r2_iso), fill = "#858585", colour = "#858585", alpha = 0.5, linetype = "dashed") +
  geom_density(data = as.data.frame(r2_hyp), aes(x = r2_hyp), fill = "#348781", colour = "#348781", alpha = 0.5) +
  geom_vline(xintercept = mean(r2_iso), color = "black", linetype = "dashed") +
  geom_vline(xintercept = mean(r2_hyp), color = "black", linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Posterior"~italic("R")^"2"*"")) + 
  mytheme

cowplot::plot_grid(green_hyper_plot, obs_pred_plot, ncol = 2, align = "h", axis = "bt")


## TEMPORAL VARIATION ##--------------------------------------------------------------------------
## changes in mass per year ##
nest_sum_grouped <- nest_dat %>%
  dplyr::mutate(repro_energy = mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * mass_kg ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001) %>%
  dplyr::group_by(year, recruit) %>%
  dplyr::summarise(mass_count   = length(mass_kg),
                   mass_mean    = mean(mass_kg),
                   mass_sd      = sd(mass_kg),
                   mass_se      = mass_sd / sqrt(mass_count),
                   # mean energy
                   energy_mean  = mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * mass_mean ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001,
                   energy_mean  = mean(repro_energy),
                   # mean clutch
                   clutch_count = length(sum_eggs),
                   clutch_mean  = mean(sum_eggs),
                   clutch_sd    = sd(sum_eggs),
                   clutch_se    = clutch_sd / sqrt(clutch_count),
                   # total
                   clutch_sum   = sum(sum_eggs),
                   energy_sum   = sum(repro_energy))

nest_sum <- nest_dat %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mass_count   = length(mass_kg),
                   mass_mean    = mean(mass_kg),
                   mass_sd      = sd(mass_kg),
                   mass_se      = mass_sd / sqrt(mass_count),
                   hyper_output = (mean(c_mydas_slope$fecundity50_int * c_mydas_slope$eggenergy_int) * mass_mean ^ mean(c_mydas_slope$fecundity50_slope + c_mydas_slope$eggenergy_slope) * 0.001),
                   clutch_count = length(sum_eggs),
                   clutch_mean  = mean(sum_eggs),
                   clutch_sd    = sd(sum_eggs),
                   clutch_se    = clutch_sd / sqrt(clutch_count),
                   clutch_sum   = sum(sum_eggs))

# Mass-year model
massyear_model <- brm(bf(mass_kg ~ s(year, by = recruit) + recruit + (1 | ID)), 
                         data    = nest_dat,  
                         family  = gaussian(),
                         #prior   = meta_priors,
                         chains  = 4, 
                         cores   = 4, 
                         iter    = 1e4, 
                         warmup  = 5e3, 
                         control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(massyear_model)

massyear_slope <- as.data.frame(coef(massyear_model)$recruit[, , 'year'])
massyear_slope$recruit <- rownames(massyear_slope) 

# posterior check
massyear_y    <- nest_dat$mass
massyear_yrep <- posterior_predict(massyear_model, draws = 500)
predmassyear  <- bayesplot::ppc_scatter_avg(massyear_y, massyear_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(55, 165), y = c(55, 165))

eggyear_model <- brms::brm(bf(sum_eggs ~ s(year, by = recruit) + recruit + (1 + year | recruit)), 
                           data    = nest_dat,
                           #prior   = meta_priors,
                           chains  = 4, 
                           cores   = 4, 
                           iter    = 1e4, 
                           warmup  = 5e3,
                           control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(eggyear_model)

eggyear_slope <- as.data.frame(coef(eggyear_model)$recruit[, , 'year'])
eggyear_slope$recruit <- rownames(massyear_slope) 

# posterior check
eggyear_y      <- nest_dat$sum_eggs
eggyear_yrep   <- posterior_predict(eggyear_model, draws = 500)
predeggyear    <- bayesplot::ppc_scatter_avg(eggyear_y, eggyear_yrep, alpha = 0.7) + geom_abline(slope = 1, linetype = "dashed") + lims(x = c(0, 1500), y = c(0, 1500))

cowplot::plot_grid(predmassyear, predeggyear, ncol = 2, align = "hv", axis = "tblr", labels = c("a", "b"))

## Construct Figure 4 ##
# Fig 4a - Count
count_plot <- ggplot() +
  geom_line(data = nest_sum, aes(x = year, y = mass_count)) +
  geom_point(data = nest_sum, aes(x = year, y = mass_count), size = 3) +
  ylab("Number of females") + xlab("Year") +
  mytheme

count_inset <- ggplot() +
  geom_line(data = nest_sum_grouped, aes(x = year, y = mass_count, colour = recruit, linetype = recruit)) +
  geom_point(data = nest_sum_grouped, aes(x = year, y = mass_count, colour = recruit, shape = recruit), size = 3) + 
  ylab("Number of females") + xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  mytheme

# Fig 4b - Fecundity count
totalclutch_plot <- ggplot() +
  geom_line(data = nest_sum, aes(x = year, y = clutch_sum)) +
  geom_point(data = nest_sum, aes(x = year, y = clutch_sum), size = 3) +
  ylab("Population Fecundity") + xlab("Year") +
  mytheme

totalclutch_inset <- ggplot() +
  geom_line(data = nest_sum_grouped, aes(x = year, y = clutch_sum, colour = recruit, linetype = recruit)) +
  geom_point(data = nest_sum_grouped, aes(x = year, y = clutch_sum, colour = recruit, shape = recruit), size = 3) + 
  ylab("Population Fecundity") + xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  mytheme

# Fig 4c - Hyper vs iso mass-scaling
hyp_iso_pred <- as.data.frame(seq(min(nest_dat$mass_kg), max(nest_dat$mass_kg))) %>% # in kg
  dplyr::rename(mass_pred = "seq(min(nest_dat$mass_kg), max(nest_dat$mass_kg))") %>%
  dplyr::mutate(hyper_50 = (mean(c_mydas_slope$fecundity50_int) * mass_pred ^ mean(c_mydas_slope$fecundity50_slope)),
                iso_50   = (mean(c_mydas_slope$fecundity50iso_int) + mean(c_mydas_slope$fecundity50iso_slope) * mass_pred),
                diff_50  = hyper_50 - iso_50)


hyp_iso_mass <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = iso_50, ymax = hyper_50), fill = "#348781", alpha = 0.2) +
  geom_point(data = nest_dat %>% filter(count_year != 1), aes(x = mass_kg, y = sum_eggs), alpha = 0.1) +
  geom_line(aes(y = iso_50), size = 1.5, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = hyper_50), size = 1.5, colour = "#348781") +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Individual Fecundity")

clutch_diff_inset <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Difference in fecundity")

# Fig 4d
hyp_iso_year <- ggplot(pred_merged, aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = pred_iso_50, ymax = pred_hyp_50), fill = "#348781", alpha = 0.2) +
  geom_point(aes(y = fecundity_obs), size = 3) +
  geom_line(aes(y = pred_iso_50), size = 1, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = pred_hyp_50), size = 1, colour = "#348781") +
  xlab("Year") +
  ylab("Population Fecundity")

year_diff_inset <- pred_merged %>%
  dplyr::mutate(diff = pred_hyp_50 - pred_iso_50) %>%
  ggplot(aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Year") +
  ylab("Difference in fecundity")

# Fig 4e - Population reproductive output
energysum_plot <- ggplot(nest_sum_grouped, aes(x = year, y = energy_sum, colour = recruit)) +
  mytheme + geom_line(aes(linetype = recruit)) +
  geom_point(aes(shape = recruit), size = 3) + 
  ylab(expression("Population "*italic(R)[i]*", (MJ)")) + 
  xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  mytheme

# Fig 4f - individual reproductive output
energymean_plot <- ggplot(nest_sum_grouped, aes(x = year, y = energy_mean, colour = recruit)) +
  geom_line(aes(linetype = recruit)) +
  geom_point(aes(shape = recruit), size = 3) + 
  ylab(expression("Individual "*italic(R)[i]*", (MJ)")) +
  xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  mytheme


# group plots
cowplot::plot_grid(count_plot, totalclutch_plot, hyp_iso_mass, hyp_iso_year, energysum_plot + theme(legend.position = "none"), energymean_plot + theme(legend.position = "none"),
          labels = c('a', 'b', 'c', 'd', 'e', 'f'), ncol = 2, align = "v", axis = 'lr') 
  

## SUPPLEMENTARY FIGURES ## -----------------------------------------------------
## Plot adult CCL and clutch frequency and mean clutch size ##
clutch_feq_iso_model <- brm(count_year ~ mass_kg, data = nest_dat, family  = gaussian(), chains = 4, cores = 4, iter = 5e3, warmup  = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 20))
clutch_feq_hyp_model <- brm(log(count_year) ~ lnMass, data = nest_dat, family  = gaussian(), chains = 4, cores = 4, iter = 5e3, warmup  = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 20))

bayes_R2(clutch_feq_iso_model)
bayes_R2(clutch_feq_hyp_model)

clutch_feq_iso_eff <- ggeffects::ggpredict(clutch_feq_iso_model, terms = "mass_kg")
clutch_feq_hyp_eff <- ggeffects::ggpredict(clutch_feq_hyp_model, terms = "lnMass")

clutch_freq_plot <- ggplot(nest_dat, aes(x = mass_kg, y = count_year)) +
  geom_point() +
  geom_line(data = clutch_feq_iso_eff, aes(x = x, y = conf.low), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_line(data = clutch_feq_iso_eff, aes(x = x, y = conf.high), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_line(data = clutch_feq_iso_eff, aes(x = x, y = predicted), linetype = "dashed", colour = "grey", size = 1) +
  geom_line(data = clutch_feq_hyp_eff, aes(x = exp(x), y = conf.low), colour = "#348781", size = 0.5) + 
  geom_line(data = clutch_feq_hyp_eff, aes(x = exp(x), y = conf.high), colour = "#348781", size = 0.5) + 
  geom_line(data = clutch_feq_hyp_eff, aes(x = exp(x), y = predicted), colour = "#348781", size = 1) +
  xlim(50, 160) +
  mytheme +
  ylab("Clutch frequency per year") +
  xlab("Female body mass (kg)")

clutch_iso_model <- brm(mean_egg ~ mass_kg, data = nest_dat, family  = gaussian(), chains = 4, cores = 4, iter = 5e3, warmup  = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 20))
clutch_hyp_model <- brm(log(mean_egg) ~ lnMass, data = nest_dat, family  = gaussian(), chains = 4, cores = 4, iter = 5e3, warmup  = 2.5e3, control = list(adapt_delta = 0.999, max_treedepth = 20))

bayes_R2(clutch_iso_model)
bayes_R2(clutch_hyp_model)

clutch_iso_eff <- ggeffects::ggpredict(clutch_iso_model, terms = "mass_kg")
clutch_hyp_eff <- ggeffects::ggpredict(clutch_hyp_model, terms = "lnMass")

clutch_plot <- ggplot(nest_dat, aes(x = mass_kg, y = mean_egg)) +
  geom_point() +
  geom_line(data = clutch_iso_eff, aes(x = x, y = conf.low), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_line(data = clutch_iso_eff, aes(x = x, y = conf.high), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_line(data = clutch_iso_eff, aes(x = x, y = predicted), linetype = "dashed", colour = "grey", size = 1) +
  geom_line(data = clutch_hyp_eff, aes(x = exp(x), y = conf.low), colour = "#348781", size = 0.5) + 
  geom_line(data = clutch_hyp_eff, aes(x = exp(x), y = conf.high), colour = "#348781", size = 0.5) + 
  geom_line(data = clutch_hyp_eff, aes(x = exp(x), y = predicted), colour = "#348781", size = 1) +
  xlim(50, 160) +
  mytheme +
  ylab("Average clutch size") +
  xlab("Female body mass (kg)")

cowplot::plot_grid(clutch_freq_plot, clutch_plot, ncol = 2, align = "hv", axis = "tblr", labels = c("a", "b"))

## Plot Mass-Year and Egg-Year modls
massyear_marg_eff <- ggeffects::ggpredict(massyear_model, terms = c("year", "recruit"))
eggyear_marg_eff  <- ggeffects::ggpredict(eggyear_model, terms = c("year", "recruit"))

massyear_plot <- ggplot() + mytheme +
  geom_ribbon(data = massyear_marg_eff, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1) +
  geom_line(data = massyear_marg_eff, aes(x = x, y = conf.low, colour = group, linetype = group), alpha = 0.5, size = 0.5) + 
  geom_line(data = massyear_marg_eff, aes(x = x, y = conf.high, colour = group, linetype = group), alpha = 0.5, size = 0.5) + 
  geom_line(data = massyear_marg_eff, aes(x = x, y = predicted, colour = group, linetype = group), size = 1) +
  geom_errorbar(data = nest_sum_grouped, aes(x = year, ymin = mass_mean-mass_se, ymax = mass_mean+mass_se, colour = recruit), width = 0.2, size = 1) +
  geom_point(data = nest_sum_grouped, aes(x = year, y = mass_mean, colour = recruit, shape = recruit), size = 3) + 
  xlim(1992, 2020) +
  ylab("Female mass, Mi, (kg)") + xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  scale_fill_manual(values = c("#9AD4A8","#14505C"))

eggyear_plot <- ggplot() + mytheme +
  geom_ribbon(data = eggyear_marg_eff, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1) +
  geom_line(data = eggyear_marg_eff, aes(x = x, y = conf.low, colour = group, linetype = group), alpha = 0.5, size = 0.5) + 
  geom_line(data = eggyear_marg_eff, aes(x = x, y = conf.high, colour = group, linetype = group), alpha = 0.5, size = 0.5) + 
  geom_line(data = eggyear_marg_eff, aes(x = x, y = predicted, colour = group, linetype = group), size = 1) +
  geom_errorbar(data = nest_sum_grouped, aes(x = year, ymin = clutch_mean-clutch_se, ymax = clutch_mean+clutch_se, colour = recruit), width = 0.2, size = 1) +
  geom_point(data = nest_sum_grouped, aes(x = year, y = clutch_mean, colour = recruit, shape = recruit), size = 3) + 
  xlim(1992, 2020) +
  ylab("Individual fecundity (total clutch size)") + xlab("Year") +
  scale_colour_manual(values = c("#9AD4A8","#14505C")) +
  scale_fill_manual(values = c("#9AD4A8","#14505C"))

cowplot::plot_grid(massyear_plot, eggyear_plot, ncol = 2, align = "hv", axis = "tblr", labels = c("a", "b"))

## Hyper vs Iso relationships Quantiles ##
head(c_mydas_slope)

hyp_iso_pred <- hyp_iso_pred %>%
  dplyr::mutate(iso_10   = (mean(c_mydas_slope$fecundity10iso_int) + mean(c_mydas_slope$fecundity10iso_slope) * mass_pred),
                iso_90   = (mean(c_mydas_slope$fecundity90iso_int) + mean(c_mydas_slope$fecundity90iso_slope) * mass_pred),
                hyper_10 = (mean(c_mydas_slope$fecundity10_int) * mass_pred ^ mean(c_mydas_slope$fecundity10_slope)),
                hyper_90 = (mean(c_mydas_slope$fecundity90_int) * mass_pred ^ mean(c_mydas_slope$fecundity90_slope)),
                diff_10  = hyper_10 - iso_10,
                diff_90  = hyper_90 - iso_90)

# Hyper vs iso mass-scaling
hyp_iso_10 <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = iso_10, ymax = hyper_10), fill = "#348781", alpha = 0.2) +
  geom_point(data = nest_dat %>% filter(count_year != 1), aes(x = mass_kg, y = sum_eggs), alpha = 0.1) +
  geom_line(aes(y = iso_10), size = 1.5, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = hyper_10), size = 1.5, colour = "#348781") +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Individual Fecundity (10th quantile)")

clutch_10_inset <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff_10), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff_10), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Difference in fecundity")

hyp_iso_90 <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = iso_90, ymax = hyper_90), fill = "#348781", alpha = 0.2) +
  geom_point(data = nest_dat %>% filter(count_year != 1), aes(x = mass_kg, y = sum_eggs), alpha = 0.1) +
  geom_line(aes(y = iso_90), size = 1.5, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = hyper_90), size = 1.5, colour = "#348781") +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Individual Fecundity (90th quantile)")

clutch_90_inset <- ggplot(hyp_iso_pred, aes(x = mass_pred)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff_90), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff_90), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab(expression("Female body mass, "*italic(M)[i]*", (kg)")) +
  ylab("Difference in fecundity")

hyp_iso_year_10 <- ggplot(pred_merged, aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = pred_iso_10, ymax = pred_hyp_10), fill = "#348781", alpha = 0.2) +
  geom_point(aes(y = fecundity_obs), size = 3) +
  geom_line(aes(y = pred_iso_10), size = 1, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = pred_hyp_10), size = 1, colour = "#348781") +
  xlab("Year") +
  ylab("Population Fecundity (10th quantile)")

year_10_inset <- pred_merged %>%
  dplyr::mutate(diff = pred_hyp_10 - pred_iso_10) %>%
  ggplot(aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Year") +
  ylab("Difference in fecundity")

hyp_iso_year_90 <- ggplot(pred_merged, aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = pred_iso_90, ymax = pred_hyp_90), fill = "#348781", alpha = 0.2) +
  geom_point(aes(y = fecundity_obs), size = 3) +
  geom_line(aes(y = pred_iso_90), size = 1, colour = "grey", linetype = "dashed") +
  geom_line(aes(y = pred_hyp_90), size = 1, colour = "#348781") +
  xlab("Year") +
  ylab("Population Fecundity (90th quantile)")

year_90_inset <- pred_merged %>%
  dplyr::mutate(diff = pred_hyp_90 - pred_iso_90) %>%
  ggplot(aes(x = year)) + mytheme +
  geom_ribbon(aes(ymin = 0, ymax = diff), fill = "#348781", alpha = 0.2) +
  geom_line(aes(y = diff), size = 1, colour = "#348781") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Year") +
  ylab("Difference in fecundity")

cowplot::plot_grid(hyp_iso_10, hyp_iso_year_10, hyp_iso_90, hyp_iso_year_90,
                   labels = c('a', 'b', 'c', 'd'), ncol = 2, align = "v", axis = 'lr') 
