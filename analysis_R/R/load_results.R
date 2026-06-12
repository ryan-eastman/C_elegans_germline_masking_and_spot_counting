# Load germquant tidy outputs from a (possibly mirrored NAS) results tree.
# Every per-image CSV already carries the shared metadata block, so a plain row-bind
# yields analysis-ready long tables.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(purrr); library(stringr)
})

# read + row-bind every "*__<suffix>.csv" under results_dir
read_tidy <- function(results_dir, suffix) {
  files <- list.files(results_dir, pattern = paste0("__", suffix, "\\.csv$"),
                      recursive = TRUE, full.names = TRUE)
  files <- files[!str_detect(files, "/_done/")]
  if (length(files) == 0) return(tibble())
  map_dfr(files, ~ suppressMessages(read_csv(.x, show_col_types = FALSE)))
}

load_results <- function(results_dir) {
  out <- list(
    nuclei         = read_tidy(results_dir, "nuclei"),
    sc_per_nucleus = read_tidy(results_dir, "sc_per_nucleus"),
    sc_tracks      = read_tidy(results_dir, "sc_tracks"),
    foci           = read_tidy(results_dir, "foci"),
    zones          = read_tidy(results_dir, "zones"),
    image_summary  = read_tidy(results_dir, "image_summary")
  )
  # canonical factor ordering for plots
  relevel_design <- function(df) {
    if (nrow(df) == 0) return(df)
    if ("treatment" %in% names(df))
      df$treatment <- factor(df$treatment, levels = c("control", "heat"))
    if ("germ_cell" %in% names(df))
      df$germ_cell <- factor(df$germ_cell, levels = c("oocyte", "spermatocyte"))
    df
  }
  lapply(out, relevel_design)
}

# join SC-per-nucleus fragmentation onto the nuclei table when present
nuclei_with_sc <- function(res) {
  n <- res$nuclei
  if (nrow(n) == 0 || nrow(res$sc_per_nucleus) == 0) return(n)
  sc <- res$sc_per_nucleus |>
    select(image_id, nucleus_id, n_fragments, sc_total_length_um,
           sc_mean_fragment_um, sc_longest_fragment_um)
  left_join(n, sc, by = c("image_id", "nucleus_id"), suffix = c("", "_sc"))
}
