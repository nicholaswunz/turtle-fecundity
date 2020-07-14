my_theme <- function() {
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.8),
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent",
                                             colour = NA))
}

wrap_grob <- function(...) {
  grid::grobTree(grid::textGrob(...))
}

make_distributions <- function(dest, fig_out_folder, ...) {
  ggplot2::ggsave(dest, distribution_figure(...), device = "pdf",
                  width = 11, height = 9, units = "in",
                  onefile = FALSE)
}

distribution_figure <- function(che_myd_data, rep_out_post) {
  che_myd_data <- che_myd_data %>%
    dplyr::mutate(year = as.factor(year))

  rep_out <- rep_out_post %>%
    dplyr::mutate(year = as.factor(year),
                  quantile = as.factor(quantile),
                  rep_out = rep_out * 1e-6)

  txta <- wrap_grob("a",
                    x = grid::unit(-0.05, "npc"),
                    y = grid::unit(1.02, "npc"),
                    gp = grid::gpar(fontsize = 18,
                                    fontface = "bold"))
  txtb <- wrap_grob("b",
                    x = grid::unit(-0.05, "npc"),
                    y = grid::unit(1.02, "npc"),
                    gp = grid::gpar(fontsize = 18,
                                    fontface = "bold"))

  a <- ggplot(data = che_myd_data,
              mapping = aes(x = mass, y = year, fill = recruit)) +
    my_theme() +
    ggridges::geom_density_ridges(scale = 1,
                                  rel_min_height = 0.01,
                                  alpha = 0.5) +
    xlab("Female mass (kg)") +
    ylab("Year") +
    colorspace::scale_fill_discrete_sequential(palette = "BluGrn") +
    coord_cartesian(clip = "off") +
    annotation_custom(txta)

  b <- ggplot(data = rep_out,
              mapping = aes(x = rep_out, y = year, fill = quantile)) +
    my_theme() +
    ggridges::geom_density_ridges(scale = 1,
                                  rel_min_height = 0.01,
                                  alpha = 0.5) +
    xlab("Population reproductive output (MJ)") +
    ylab("Year") +
    colorspace::scale_fill_discrete_sequential(palette = "Heat") +
    scale_x_continuous(trans = "log10") +
    coord_cartesian(clip = "off") +
    annotation_custom(txtb)
  p <- gridExtra::arrangeGrob(a, b, ncol = 2)
  ggpubr::annotate_figure(p, top = ggpubr::text_grob(""))
}

make_r2s_hist <- function(dest, fig_out_folder, ...) {
  ggplot2::ggsave(dest, r2s_hist(...), device = "pdf",
                  width = 6, height = 6, units = "in",
                  onefile = FALSE)
}

r2s_hist <- function(plot_data) {
  r2_iso <- plot_data$r2_vals[plot_data$type == "iso"]
  r2_hyp <- plot_data$r2_vals[plot_data$type == "hyp"]
  ggplot(plot_data, aes(x = r2_vals)) +
    geom_density(data = plot_data, adjust = 2, trim = TRUE,
                 mapping = aes(x = r2_vals,
                               y = ..scaled..,
                               group = type,
                               colour = type,
                               linetype = type,
                               fill = type),
                 alpha = 0.5,
                 lwd = 0.8,
                 show.legend = FALSE) +
    geom_vline(xintercept = median(r2_iso), linetype = 2,
               colour = "grey30", size = 0.3) +
    geom_vline(xintercept = median(r2_hyp), linetype = 2,
               colour = "grey30", size = 0.3) +
    scale_colour_manual(values = c("tomato", "dodgerblue2")) +
    scale_fill_manual(values = c("tomato", "dodgerblue2")) +
    theme_classic() +
    labs(y = "Scaled [0,1] posterior density",
         x = substitute("Posterior " * italic("R"^2))) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0))
}
