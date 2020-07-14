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
