#!/usr/bin/env Rscript
# =============================================================================
# analyze_coloc.R — SYP<->PGL-1 perinuclear-shell colocalization: grouped
# comparison (sex x treatment x genotype) + agreement with Imaris. ggplot2.
#
# Reads the pipeline's per-gonad `*__coloc.csv` files, keeps the HEADLINE
# `shell_voxel` row (SYP<->PGL-1 voxel coloc in the lamin-defined perinuclear
# shell), parses metadata from the image id, and makes the figures + stats.
# Open in RStudio/Positron and Source, or: Rscript scripts/analyze_coloc.R
# =============================================================================

suppressPackageStartupMessages({
  ok <- require(tidyverse, quietly = TRUE)
})
if (!ok) stop("Install tidyverse:  install.packages('tidyverse')")

# ---- CONFIG: edit these paths --------------------------------------------------
result_dirs <- c(
  "C:/Users/ryane/ccw77_fullres",  # ccw77 20260708 (10 gonads)
  "C:/Users/ryane/ccw77_0622",     # ccw77 20260622 (+ herm replicates)
  "C:/Users/ryane/n2_0625"         # N2 wild-type control
)
imaris_csv <- "C:/Users/ryane/coloc_analysis/imaris_gt.csv"  # image_id,imaris_pearson (see below)
out_dir    <- "C:/Users/ryane/coloc_analysis"
# --------------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- read the shell_voxel coloc row from every gonad ----
coloc_files <- result_dirs |>
  keep(dir.exists) |>
  map(\(d) list.files(d, "__coloc\\.csv$", full.names = TRUE, recursive = TRUE)) |>
  list_c()
if (length(coloc_files) == 0) stop("No __coloc.csv found under result_dirs — is the batch done?")

read_shell <- function(f) {
  df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  df <- dplyr::filter(df, sc_operand == "shell_voxel")
  if (nrow(df) == 0) return(NULL)          # gonads run before the shell metric existed
  # keep only what we use (metadata is parsed from image_id below) and force stable types,
  # so binding across files never clashes on an unused column (e.g. replicate int vs chr).
  dplyr::transmute(df,
    image_id = as.character(image_id), region = as.character(region),
    pearson_r = as.numeric(pearson_r),
    manders_m1 = as.numeric(manders_m1), manders_m2 = as.numeric(manders_m2),
    source_file = basename(f))
}
raw <- purrr::map(coloc_files, read_shell) |> purrr::list_rbind()

# ---- parse sex / treatment / genotype / date from image_id (robust; does not
#      depend on the pipeline's own filename parsing) ----
dat <- raw |>
  mutate(
    date      = str_extract(image_id, "^\\d{6,8}"),
    genotype  = str_to_lower(str_extract(image_id, regex("ccw77|n2", ignore_case = TRUE))),
    treatment = if_else(str_detect(image_id, regex("nohs", ignore_case = TRUE)), "noHS", "HS"),
    sex       = if_else(str_detect(image_id, regex("herm", ignore_case = TRUE)), "herm", "male"),
    rep       = str_extract(image_id, "\\d+$"),
    treatment = factor(treatment, c("noHS", "HS")),
    sex       = factor(sex, c("herm", "male")),
    genotype  = factor(genotype, c("n2", "ccw77"))
  ) |>
  transmute(image_id, date, genotype, treatment, sex, rep,
            shell_pearson = pearson_r, shell_m1 = manders_m1, shell_m2 = manders_m2,
            region, source_file) |>
  arrange(genotype, treatment, sex, rep)

readr::write_csv(dat, file.path(out_dir, "coloc_combined.csv"))
message(sprintf("Loaded %d gonads across %d experiment date(s).",
                nrow(dat), dplyr::n_distinct(dat$date)))

# ---- Fig 1: shell Pearson by sex, faceted genotype x treatment ----
p1 <- ggplot(dat, aes(sex, shell_pearson, fill = sex)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.35) +
  geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.85) +
  facet_grid(genotype ~ treatment, labeller = label_both) +
  scale_fill_manual(values = c(herm = "#E69F00", male = "#0072B2")) +
  labs(title = "SYP↔PGL-1 perinuclear colocalization",
       subtitle = "shell Pearson (lamin-defined perinuclear cytoplasm; SC ribbon excluded)",
       x = NULL, y = "shell Pearson") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out_dir, "fig1_shell_pearson_by_group.png"), p1,
       width = 7, height = 6, dpi = 200)

# ---- stats: male vs herm within each genotype x treatment ----
ccc <- function(x, y) {                       # Lin's concordance correlation coefficient
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}
stat_tab <- dat |>
  group_by(genotype, treatment) |>
  summarise(
    n_herm = sum(sex == "herm"), n_male = sum(sex == "male"),
    herm_mean = mean(shell_pearson[sex == "herm"]),
    male_mean = mean(shell_pearson[sex == "male"]),
    diff_male_minus_herm = male_mean - herm_mean,
    p_wilcox = tryCatch(
      if (n_herm > 0 && n_male > 0)
        suppressWarnings(wilcox.test(shell_pearson ~ sex)$p.value) else NA_real_,
      error = \(e) NA_real_),
    .groups = "drop")
readr::write_csv(stat_tab, file.path(out_dir, "stats_male_vs_herm.csv"))
message("male-vs-herm summary:"); print(stat_tab)

# ---- Fig 2: agreement with Imaris (only if the GT csv exists) ----
# imaris_gt.csv needs columns: image_id, imaris_pearson  (fill it as you do Imaris).
# If it doesn't exist yet, we write a template listing every gonad for you to fill in.
if (!file.exists(imaris_csv)) {
  readr::write_csv(tibble(image_id = dat$image_id, imaris_pearson = NA_real_), imaris_csv)
  message("Wrote a template to fill: ", imaris_csv)
} else {
  im <- suppressMessages(readr::read_csv(imaris_csv, show_col_types = FALSE)) |>
    filter(!is.na(imaris_pearson))
  ag <- inner_join(dat, im, by = "image_id")
  if (nrow(ag) >= 2) {
    cval <- ccc(ag$shell_pearson, ag$imaris_pearson)
    lims <- range(c(ag$shell_pearson, ag$imaris_pearson), na.rm = TRUE)
    p2 <- ggplot(ag, aes(imaris_pearson, shell_pearson)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(aes(colour = sex, shape = genotype), size = 3.2) +
      scale_colour_manual(values = c(herm = "#E69F00", male = "#0072B2")) +
      coord_equal(xlim = lims, ylim = lims) +
      labs(title = sprintf("Pipeline vs Imaris (CCC = %.2f, n = %d gonads)", cval, nrow(ag)),
           x = "Imaris Pearson", y = "pipeline shell Pearson") +
      theme_bw(base_size = 12)
    ggsave(file.path(out_dir, "fig2_imaris_agreement.png"), p2,
           width = 5.5, height = 5, dpi = 200)
    message(sprintf("Imaris agreement: CCC = %.3f over %d gonads.", cval, nrow(ag)))
  } else {
    message("imaris_gt.csv has < 2 filled rows — add more imaris_pearson values.")
  }
}

message("Done. Outputs in: ", out_dir)
