#!/usr/bin/env Rscript
# Build all figures + a summary table from a germquant results folder.
#   Rscript analysis_R/analyze.R <results_dir> [out_dir]
# Or open germquant_analysis.Rproj in Positron/RStudio and source this interactively.
suppressPackageStartupMessages({ library(ggplot2); library(readr) })

here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ...)
src_dir <- tryCatch(here("R"), error = function(e) "analysis_R/R")
if (!dir.exists(src_dir)) src_dir <- "analysis_R/R"
source(file.path(src_dir, "load_results.R"))
source(file.path(src_dir, "figures.R"))

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
out_dir     <- if (length(args) >= 2) args[2] else file.path(results_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading results from: ", results_dir)
res <- load_results(results_dir)
nsc <- nuclei_with_sc(res)
message(sprintf("  %d nuclei across %d images", nrow(res$nuclei),
                length(unique(res$nuclei$image_id))))

save_plot <- function(p, name, w = 6, h = 4) {
  if (is.null(p)) return(invisible())
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 200)
  message("  wrote ", name)
}

save_plot(plot_sc_fragmentation(nsc), "sc_fragmentation.png")
save_plot(plot_sc_length(nsc),        "sc_length.png")
save_plot(plot_rad51(res$nuclei),     "rad51_foci.png")
save_plot(plot_zone_length(res$zones), "zone_length.png", w = 6, h = 5)

st <- summary_table(nsc)
if (nrow(st) > 0) {
  write_csv(st, file.path(out_dir, "group_summary.csv"))
  message("\nGroup summary (pachytene nuclei):"); print(st)
}
message("\nDone -> ", out_dir)
