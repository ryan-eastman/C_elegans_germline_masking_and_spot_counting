# Publication-style figures for the heat × sex SC/foci/zone readouts.
# x = treatment (control/heat), facet = germ_cell (oocyte/spermatocyte), violin+box+jitter.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

GQ_THEME <- theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none")
GQ_FILL <- scale_fill_manual(values = c(control = "#4C72B0", heat = "#C44E52"))

# generic grouped violin for one per-nucleus metric
.violin <- function(df, y, ylab, title, pachytene_only = TRUE) {
  if (!all(c(y, "treatment", "germ_cell") %in% names(df))) {
    message("  [skip] '", title, "' — missing column ", y); return(NULL)
  }
  if (pachytene_only && "zone_call" %in% names(df))
    df <- dplyr::filter(df, zone_call == "pachytene")
  df <- dplyr::filter(df, !is.na(.data[[y]]), !is.na(treatment), !is.na(germ_cell))
  if (nrow(df) == 0) { message("  [skip] '", title, "' — no rows"); return(NULL) }
  ggplot(df, aes(treatment, .data[[y]], fill = treatment)) +
    geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.08, size = 0.3, alpha = 0.25) +
    facet_wrap(~ germ_cell) +
    labs(x = NULL, y = ylab, title = title) +
    GQ_FILL + GQ_THEME
}

# THE heat phenotype: SC fragmentation (fragments / nucleus)
plot_sc_fragmentation <- function(nuclei_sc)
  .violin(nuclei_sc, "n_fragments", "SC fragments / nucleus",
          "SC fragmentation (heat fragments the SC in spermatocytes)")

plot_sc_length <- function(nuclei_sc)
  .violin(nuclei_sc, "sc_total_length_um", "Total SC length / nucleus (µm)",
          "Synaptonemal complex length")

plot_rad51 <- function(nuclei)
  .violin(nuclei, "n_foci", "RAD-51 foci / nucleus",
          "RAD-51 DNA-damage foci")

# zone lengths (one row per zone per germline) — not per-nucleus
plot_zone_length <- function(zones) {
  if (nrow(zones) == 0 || !all(c("zone", "length_um", "treatment", "germ_cell") %in% names(zones))) {
    message("  [skip] zone length — no/empty zones table"); return(NULL)
  }
  ggplot(zones, aes(treatment, length_um, fill = treatment)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.6) +
    facet_grid(zone ~ germ_cell, scales = "free_y") +
    labs(x = NULL, y = "Zone length (µm)", title = "Transition-zone & pachytene length") +
    GQ_FILL + GQ_THEME
}

# group means ± sd — the numbers behind the plots
summary_table <- function(nuclei_sc) {
  metrics <- intersect(c("n_fragments", "sc_total_length_um", "n_foci"), names(nuclei_sc))
  if (length(metrics) == 0) return(tibble())
  df <- nuclei_sc
  if ("zone_call" %in% names(df)) df <- dplyr::filter(df, zone_call == "pachytene")
  df |>
    group_by(germ_cell, treatment) |>
    summarise(n_nuclei = dplyr::n(),
              across(all_of(metrics), list(mean = ~mean(.x, na.rm = TRUE),
                                           sd = ~sd(.x, na.rm = TRUE))),
              .groups = "drop")
}
