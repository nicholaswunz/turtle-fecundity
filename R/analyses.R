read_file <- function(filename, ...) {
  read.csv(filename, header = TRUE, stringsAsFactors = FALSE, ...)
}

clean_turtle_data <- function(...) {
  read_file(...) %>%
    dplyr::filter(species != "") %>%
    dplyr::mutate(ln_mass_g = log(pred.mass),
                  ln_clutch_size = log(pred.clutch),
                  ln_egg_vol = log(egg.volume),
                  animal = as.factor(species))
}

clean_energy_data <- function(...) {
  read_file(...) %>%
    dplyr::filter(species == "Chelonia_mydas") %>%
    dplyr::mutate(ln_energy_j = log(egg.energy),
                  ln_egg_mass_g = log(egg.mass))
}

clean_che_myd_data <- function(...) {
  read_file(...) %>%
    dplyr::mutate(ln_clutch_size = log(sum_n.of.eggs),
                  ln_mass_g = log(mass * 1e3))
}

run_phylo_model <- function(data, y_var, phylo_cor) {
  data$y <- data[[y_var]]
  brms:::brm(y ~ ln_mass_g + (1 + ln_mass_g | species) +
               (1 | gr(animal, cov = phylo_cor)),
             data = data, family = gaussian(),
             data2 = list(phylo_cor = phylo_cor),
             prior = c(brms::prior(normal(0, 10), "b"),
                       brms::prior(normal(0, 10), "Intercept"),
                       brms::prior(student_t(3, 0, 20), "sd"),
                       brms::prior(student_t(3, 0, 20), "sigma")),
             chains = 4, cores = 4, iter = 5e3,
             warmup = 2.5e3, control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))
}

run_energy_model <- function(data) {
  brms:::brm(ln_energy_j ~ ln_egg_mass_g,
             data = data, family = gaussian(),
             prior = c(brms::prior(normal(0, 10), "b"),
                       brms::prior(normal(0, 10), "Intercept"),
                       brms::prior(student_t(3, 0, 20), "sigma")),
             chains = 4, cores = 4, iter = 5e3,
             warmup = 2.5e3, control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))
}

extract_post_params <- function(model) {
  post_data <- brms::posterior_samples(model)
  species <- sort(unique(model$data$species))
  out_list <- vector(mode = "list", length = length(species))
  mean_int <- "b_Intercept"
  mean_slp <- "b_ln_mass_g"
  for (j in seq_along(species)) {
    out_list[[j]] <- data.frame(matrix(NA, nrow(post_data), 2))
    dat_names <- paste0(c("Intercept_", "Slope_"), species[j])
    names(out_list[[j]]) <- dat_names
    phylo_int <- paste0("r_animal[", species[j], ",Intercept]")
    sp_int <- paste0("r_species[", species[j], ",Intercept]")
    sp_slp <- paste0("r_species[", species[j], ",ln_mass_g]")
    for (i in seq_len(nrow(post_data))) {
      out_list[[j]][i, dat_names[1]] <- exp(post_data[i, mean_int] +
                                              post_data[i, phylo_int] +
                                              post_data[i, sp_int])
      out_list[[j]][i, dat_names[2]] <- post_data[i, mean_slp] +
                                          post_data[i, sp_slp]
    }
  }
  do.call("cbind.data.frame", out_list)
}

melted_posteriors <- function(param_data) {
  set.seed(10)
  param_data <- dplyr::sample_n(param_data, 4000)
  int_targets <- grep("Intercept", names(param_data), value = TRUE)
  slp_targets <- grep("Slope", names(param_data), value = TRUE)
  # merge fecundity intercept data
  intercept <- reshape2::melt(param_data, measure.vars = int_targets,
                              variable.name = "species",
                              value.name = "intercept")
  # merge fecundity slope data
  slope <- reshape2::melt(param_data, measure.vars = slp_targets,
                          variable.name = "species",
                          value.name = "slope")

  data.frame(species = intercept$species,
             intercept = intercept$intercept,
             slope = slope$slope,
             stringsAsFactors = FALSE)
}

extract_and_melt <- function(...) {
  extract_post_params(...) %>%
    melted_posteriors() %>%
    dplyr::mutate(species = as.factor(substr(species, 11, 100)))
}

make_tot_clutch_vol_data <- function(clutch_params, egg_params) {
  data.frame(clutch_params, egg_params) %>%
    dplyr::mutate(intercept = intercept * intercept.1,
                  slope = slope + slope.1)
}

extract_energy_params <- function(model) {
  brms::posterior_samples(model) %>%
    dplyr::sample_n(4000) %>%
    dplyr::mutate(species = "Chelonia_mydas",
                  intercept = exp(b_Intercept),
                  slope = b_ln_egg_mass_g) %>%
    dplyr::select(species, intercept, slope)
}

run_quantile_regression <- function(data, quantile) {
  brms::brm(brms::bf(ln_clutch_size ~ ln_mass_g, quantile = quantile),
            data = data, family = asym_laplace(),
            prior = c(brms::prior(normal(0, 10), "b"),
                      brms::prior(normal(0, 10), "Intercept"),
                      brms::prior(student_t(3, 0, 20), "sigma")),
             chains = 4, cores = 4, iter = 5e3,
             warmup = 2.5e3, control = list(adapt_delta = 0.999,
                                            max_treedepth = 15))
}

run_che_myd_models <- function(data) {
  list(quantile_10 = run_quantile_regression(data, 0.1),
       quantile_50 = run_quantile_regression(data, 0.5),
       quantile_90 = run_quantile_regression(data, 0.9),
       mean = brms::brm(ln_clutch_size ~ ln_mass_g,
                        data = data, family = gaussian(),
                        prior = c(brms::prior(normal(0, 10), "b"),
                                  brms::prior(normal(0, 10), "Intercept"),
                                  brms::prior(student_t(3, 0, 20), "sigma")),
                         chains = 4, cores = 4, iter = 5e3,
                         warmup = 2.5e3, control = list(adapt_delta = 0.999,
                                                        max_treedepth = 15)))
}

sub_param_data <- function(data, tag) {
  data <- data %>%
    dplyr::filter(species == "Chelonia_mydas")
  names(data)[2:3] <- paste0(tag, "_", names(data)[2:3])
  data
}

concatenate_posteriors <- function(eng_params, egg_params, all_models) {
  out <- plyr::llply(names(all_models), function(x, all_models) {
    posts <- brms::posterior_samples(all_models[[x]]) %>%
      dplyr::sample_n(4000) %>%
      dplyr::mutate(intercept = exp(b_Intercept),
                    slope = b_ln_mass_g) %>%
      dplyr::select(intercept, slope)
    names(posts) <- paste0(x, "_", names(posts))
    posts
  }, all_models = all_models) %>%
    do.call("cbind.data.frame", args = .)
  egg_params <- sub_param_data(egg_params, "eggvol")
  eng_params <- sub_param_data(eng_params,
                               "eggeng") %>%
    dplyr::select(-species)
  # assume 1 cm3 = 1 g (from egg energy model)
  cbind(data.frame(egg_params, eng_params), out) %>%
    dplyr::mutate(repro_int_10 = quantile_10_intercept *
                                   eggeng_intercept *
                                   (eggvol_intercept^eggeng_slope),
                  repro_slp_10 = quantile_10_slope +
                                   (eggvol_slope * eggeng_slope),
                  repro_int_50 = quantile_50_intercept *
                                   eggeng_intercept *
                                   (eggvol_intercept^eggeng_slope),
                  repro_slp_50 = quantile_50_slope +
                                   (eggvol_slope * eggeng_slope),
                  repro_int_90 = quantile_90_intercept *
                                   eggeng_intercept *
                                   (eggvol_intercept^eggeng_slope),
                  repro_slp_90 = quantile_90_slope +
                                   (eggvol_slope * eggeng_slope))
}

calc_rep_out_post <- function(posteriors, data) {
  yr_summary <- data.frame()
  for (i in seq_len(nrow(data))) {
    data$repro_10 <- posteriors$repro_int_10[i] *
      (data$mass * 1e3)^posteriors$repro_slp_10[i]
    data$repro_50 <- posteriors$repro_int_50[i] *
      (data$mass * 1e3)^posteriors$repro_slp_50[i]
    data$repro_90 <- posteriors$repro_int_90[i] *
      (data$mass * 1e3)^posteriors$repro_slp_90[i]
    tmp <- data %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(pop_rep_out_10 = sum(repro_10),
                       pop_rep_out_50 = sum(repro_50),
                       pop_rep_out_90 = sum(repro_90)) %>%
      dplyr::ungroup() %>%
      reshape2::melt(measure.vars = c("pop_rep_out_10",
                                      "pop_rep_out_50",
                                      "pop_rep_out_90"),
                     variable.name = "quantile",
                     value.name = "rep_out") %>%
      dplyr::mutate(quantile = as.numeric(gsub("pop_rep_out_", "0.", quantile)),
                    mcmc_draw = i)
    yr_summary <- rbind(yr_summary, tmp)
  }
  yr_summary
}
