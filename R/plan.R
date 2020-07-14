plan <- drake::drake_plan(
  # data -----------------------------------------------
  turtle_data = clean_turtle_data("data/reproductive output.csv"),
  energy_data = clean_energy_data("data/egg energy.csv"),
  che_myd_data = clean_che_myd_data("data/nesting year.csv"),

  # analysis -------------------------------------------
  species = sort(unique(turtle_data$species)),
  taxa = rotl::tnrs_match_names(names = species),
  tree = rotl::tol_induced_subtree(ott_ids = rotl::ott_id(taxa),
                                   label_format = "name"),
  phylo_tree = ape::compute.brlen(tree, method = "Grafen", power = 1),
  phylo_cor = ape::vcv(phylo_tree, cor = TRUE),
  clutch_model = run_phylo_model(turtle_data, "ln_clutch_size", phylo_cor),
  clutch_params = extract_and_melt(clutch_model),
  egg_model = run_phylo_model(turtle_data, "ln_egg_vol", phylo_cor),
  egg_params = extract_and_melt(egg_model),
  tot_clutch = make_tot_clutch_vol_data(clutch_params, egg_params),
  che_myd_energy_model = run_energy_model(energy_data),
  che_myd_energy_params = extract_energy_params(che_myd_energy_model),
  che_myd_clutch_models = run_che_myd_models(che_myd_data),
  che_myd_clutch_model_iso = run_quantile_regression_iso(che_myd_data, 0.5),
  che_myd_clutch_params = concatenate_posteriors(che_myd_energy_params,
                                                 egg_params,
                                                 che_myd_clutch_models),
  rep_out_post = calc_rep_out_post(che_myd_clutch_params, che_myd_data),
  models_r2s = compare_r2s(che_myd_data, che_myd_clutch_model_iso,
                           che_myd_clutch_models),

  # figures --------------------------------------------
  fig_out_folder = dir.create("output/figures/",
                              recursive = TRUE,
                              showWarnings = FALSE),
  fig_3d = make_r2s_hist(file_out("output/figures/fig_3d.pdf"),
                         fig_out_folder, models_r2s),
  ed_fig_4 = make_distributions(file_out("output/figures/ed_fig_4.pdf"),
                                fig_out_folder, che_myd_data, rep_out_post)
)
