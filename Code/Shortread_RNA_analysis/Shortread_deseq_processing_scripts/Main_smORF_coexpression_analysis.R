# =============================================================================
# Main-ORF / smORF co-expression analysis (coupled vs non-coupled routing)
#
# For each smORF / main-ORF pair in the MASTER annotation table, this script:
#   1. Pulls batch-corrected VST expression for the smORF and its main ORF
#      from the DESeq2 object produced by `RNA_differential_expression.R`.
#   2. Computes Pearson r in each diagnosis group (default: Control vs AD)
#      and Δr (= r_AD − r_Control).
#   3. Produces a per-pair triptych figure
#         [scatter w/ group-coloured fits | main-ORF violin | smORF violin]
#      and writes it to:
#         <OUT_ROOT>/coupled/      if |Δr| <= DELTA_R_ROUTE_THRESHOLD
#         <OUT_ROOT>/non_coupled/  otherwise
#   4. Saves a tidy summary CSV plus an RDS of the per-pair fit-line data
#      that can be reloaded to rebuild the global "landscape" figure.
#
# Generalizable: change CONFIG paths/labels to switch cohorts. Consumes the
# `<prefix>_dds.rds` and `<prefix>_vsd_corrected.rds` produced upstream.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---------------------------- CONFIG -----------------------------------------
# Paths are relative to DESEQ_DIR / DATA_DIR / OUT_ROOT. Override via env vars
# `AD_PAPER_DESEQ_DIR`, `AD_PAPER_DATA_DIR`, `AD_PAPER_COEXPR_OUT_DIR`.
DESEQ_DIR <- Sys.getenv("AD_PAPER_DESEQ_DIR",      unset = "./outputs/deseq")
DATA_DIR  <- Sys.getenv("AD_PAPER_DATA_DIR",       unset = "./data")
OUT_ROOT  <- Sys.getenv("AD_PAPER_COEXPR_OUT_DIR", unset = "./outputs/main_smorf_coexpression")

DATASET <- "ROSMAP"
PREFIX  <- tolower(DATASET)

DDS_PATH    <- file.path(DESEQ_DIR, sprintf("%s_dds.rds", PREFIX))
VSD_PATH    <- file.path(DESEQ_DIR, sprintf("%s_vsd_corrected.rds", PREFIX))
MASTER_PATH <- file.path(DATA_DIR, "Brain_Microproteins_Discovery_MASTER_dataset.csv")

GROUP_COL    <- "new_dx"
GROUP_LEVELS <- c("Control", "AD")
GROUP_LABELS <- c("CI", "SymAD")             # display labels for plots
GROUP_COLORS <- c("CI" = "black", "SymAD" = "#b70000")
COUPLED_DIR     <- file.path(OUT_ROOT, "coupled")
NON_COUPLED_DIR <- file.path(OUT_ROOT, "non_coupled")

DELTA_R_ROUTE_THRESHOLD <- 0.1
EXCLUDE_TYPES <- c("lncRNA", "psORF", "eORF")

BACKGROUND_COLORS <- c("SwissProt" = "#74a2b7", "Noncanonical" = "#ed8651")
# -----------------------------------------------------------------------------

dir.create(COUPLED_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(NON_COUPLED_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Gradient background ----------------------------
create_gradient_bg <- function(col1 = BACKGROUND_COLORS[["SwissProt"]],
                               col2 = BACKGROUND_COLORS[["Noncanonical"]],
                               col_mid = "white",
                               angle_deg = -50,
                               n = 400) {
  ang <- angle_deg * pi / 180
  dx <- cos(ang); dy <- sin(ang)
  xs <- seq(-1, 1, length.out = n); ys <- seq(-1, 1, length.out = n)
  gx <- matrix(xs, n, n, byrow = TRUE)
  gy <- matrix(rev(ys), n, n)
  proj <- gx * dx + gy * dy
  t <- (proj - min(proj)) / (max(proj) - min(proj))
  c1 <- col2rgb(col1) / 255; cm <- col2rgb(col_mid) / 255; c2 <- col2rgb(col2) / 255
  t2 <- t * 2; lo <- t2 <= 1
  r <- ifelse(lo, c1[1] + t2 * (cm[1] - c1[1]), cm[1] + (t2 - 1) * (c2[1] - cm[1]))
  g <- ifelse(lo, c1[2] + t2 * (cm[2] - c1[2]), cm[2] + (t2 - 1) * (c2[2] - cm[2]))
  b <- ifelse(lo, c1[3] + t2 * (cm[3] - c1[3]), cm[3] + (t2 - 1) * (c2[3] - cm[3]))
  grid::rasterGrob(array(c(r, g, b), dim = c(n, n, 3)),
                   width  = grid::unit(1, "npc"),
                   height = grid::unit(1, "npc"))
}

# ---------------------------- Load inputs ------------------------------------
master <- read.csv(MASTER_PATH, stringsAsFactors = FALSE)
dds      <- readRDS(DDS_PATH)
vst_expr <- readRDS(VSD_PATH)
meta_df  <- as.data.frame(colData(dds))

stopifnot(GROUP_COL %in% colnames(meta_df), length(GROUP_LEVELS) == 2)

g1 <- rownames(meta_df[meta_df[[GROUP_COL]] == GROUP_LEVELS[1], ])
g2 <- rownames(meta_df[meta_df[[GROUP_COL]] == GROUP_LEVELS[2], ])
n1 <- length(g1); n2 <- length(g2)
cat(sprintf("%s: %d, %s: %d\n", GROUP_LEVELS[1], n1, GROUP_LEVELS[2], n2))

# ---------------------------- Build pairs ------------------------------------
coord_col <- dplyr::case_when(
  "genomic_coordinate"  %in% colnames(master) ~ "genomic_coordinate",
  "genomic_coordinates" %in% colnames(master) ~ "genomic_coordinates",
  "coordinates"         %in% colnames(master) ~ "coordinates",
  TRUE ~ NA_character_
)
if (is.na(coord_col)) stop("No coordinate column in MASTER.")
if (!"RNA_Ensembl_Parent" %in% colnames(master))
  stop("RNA_Ensembl_Parent column required.")

pairs_df <- master %>%
  dplyr::mutate(
    main_id_raw  = trimws(as.character(RNA_Ensembl_Parent)),
    smorf_coord  = trimws(as.character(.data[[coord_col]])),
    smorf_symbol = trimws(as.character(gene_symbol)),
    smorf_match  = dplyr::coalesce(dplyr::na_if(smorf_coord, ""),
                                   dplyr::na_if(smorf_symbol, ""))
  ) %>%
  dplyr::filter(
    !trimws(smorf_type) %in% EXCLUDE_TYPES,
    grepl("^ENSG", main_id_raw),
    !grepl(",", main_id_raw),
    !is.na(smorf_match)
  ) %>%
  dplyr::transmute(
    smorf_id  = smorf_match,
    main_id   = main_id_raw,
    main_gene = gene_name
  ) %>%
  dplyr::distinct()

cat(sprintf("Pairs to plot: %d\n", nrow(pairs_df)))

# ---------------------------- Helpers ----------------------------------------
find_row <- function(id, expr_rownames) {
  if (id %in% expr_rownames) return(id)
  regex <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", id)
  regex <- gsub(":", ".", regex)
  m <- grep(regex, expr_rownames, value = TRUE)
  if (length(m) > 0) return(m[1]) else NA_character_
}

make_violin_panel <- function(expr_vec, bg_type, title_label) {
  diag_vec <- c(rep(GROUP_LABELS[1], n1), rep(GROUP_LABELS[2], n2))
  df_v <- data.frame(
    count     = as.numeric(expr_vec),
    diagnosis = factor(diag_vec, levels = GROUP_LABELS)
  )
  bg_col <- BACKGROUND_COLORS[[bg_type]]
  ggplot(df_v, aes(x = diagnosis, y = count, fill = diagnosis)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = scales::alpha(bg_col, 0.52), color = NA) +
    geom_violin(trim = FALSE, color = "black", linewidth = 0.6, width = 0.9) +
    geom_boxplot(width = 0.15, color = "black", fill = "white",
                 linewidth = 0.6, outlier.alpha = 0, alpha = 0.4) +
    scale_fill_manual(values = c(setNames(c("grey70", "#b70000"), GROUP_LABELS))) +
    labs(title = title_label, x = NULL, y = "VST") +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 11),
      axis.text        = element_text(size = 11, color = "black"),
      axis.title.y     = element_text(size = 12, margin = margin(r = 8)),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.45),
      legend.position  = "none"
    )
}

# ---------------------------- Main loop --------------------------------------
gradient_bg <- create_gradient_bg()
rn <- rownames(vst_expr)

results <- list(); global_lines <- list()
attempted <- nrow(pairs_df)
n_skip_missing <- n_skip_invalid <- n_coupled <- n_non_coupled <- 0

for (i in seq_len(nrow(pairs_df))) {
  smorf_raw <- pairs_df$smorf_id[i]
  main_raw  <- pairs_df$main_id[i]
  main_gene <- pairs_df$main_gene[i]

  smorf_row <- find_row(smorf_raw, rn)
  main_row  <- find_row(main_raw,  rn)
  if (is.na(smorf_row) || is.na(main_row)) {
    n_skip_missing <- n_skip_missing + 1; next
  }

  s1 <- as.numeric(vst_expr[smorf_row, g1])
  p1 <- as.numeric(vst_expr[main_row,  g1])
  s2 <- as.numeric(vst_expr[smorf_row, g2])
  p2 <- as.numeric(vst_expr[main_row,  g2])

  if (sum(stats::complete.cases(s1, p1)) < 3 ||
      sum(stats::complete.cases(s2, p2)) < 3) {
    n_skip_invalid <- n_skip_invalid + 1; next
  }

  r1 <- suppressWarnings(stats::cor(s1, p1, use = "complete.obs"))
  r2 <- suppressWarnings(stats::cor(s2, p2, use = "complete.obs"))
  if (is.na(r1) || is.na(r2)) { n_skip_invalid <- n_skip_invalid + 1; next }

  delta_r <- r2 - r1
  route   <- if (abs(delta_r) > DELTA_R_ROUTE_THRESHOLD) "non_coupled" else "coupled"

  df_plot <- data.frame(
    smorf  = c(s1, s2),
    parent = c(p1, p2),
    group  = factor(c(rep(GROUP_LABELS[1], n1), rep(GROUP_LABELS[2], n2)),
                    levels = GROUP_LABELS),
    stringsAsFactors = FALSE
  )

  global_lines[[length(global_lines) + 1]] <- data.frame(
    pair_id = sprintf("%s__%s", main_raw, smorf_raw),
    smorf   = df_plot$smorf,
    parent  = df_plot$parent,
    group   = df_plot$group,
    stringsAsFactors = FALSE
  )

  p_scatter <- ggplot(df_plot, aes(x = smorf, y = parent, color = group, fill = group)) +
    annotation_custom(gradient_bg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_point(alpha = 0.4, size = 1.4, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_color_manual(values = GROUP_COLORS) +
    scale_fill_manual(values  = GROUP_COLORS) +
    guides(fill = "none") +
    labs(
      x = sprintf("smORF: %s (VST)", smorf_raw),
      y = sprintf("Main ORF: %s (VST)", main_raw),
      title = sprintf("%s | %s", main_gene, smorf_raw),
      subtitle = sprintf("r_%s = %.3f, r_%s = %.3f, Δr = %.3f",
                         GROUP_LABELS[1], r1, GROUP_LABELS[2], r2, delta_r),
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = NA, color = NA),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text        = element_text(size = 11, color = "black"),
      axis.title       = element_text(size = 12),
      plot.title       = element_text(size = 11),
      plot.subtitle    = element_text(size = 9)
    )

  triptych <- p_scatter +
    make_violin_panel(c(p1, p2), "SwissProt",   main_gene) +
    make_violin_panel(c(s1, s2), "Noncanonical", smorf_raw) +
    patchwork::plot_layout(widths = c(1.5, 1, 1))

  out_dir  <- if (route == "non_coupled") NON_COUPLED_DIR else COUPLED_DIR
  fname    <- sprintf("%s_%s.pdf",
                      gsub("[^A-Za-z0-9_.-]", "_", main_gene),
                      gsub("[^A-Za-z0-9_.-]", "_", smorf_raw))
  ggsave(file.path(out_dir, fname), plot = triptych,
         width = 16, height = 5, device = cairo_pdf)

  if (route == "non_coupled") n_non_coupled <- n_non_coupled + 1 else n_coupled <- n_coupled + 1

  results[[length(results) + 1]] <- data.frame(
    smorf_id = smorf_raw, main_id = main_raw, main_gene_name = main_gene,
    smorf_row = smorf_row, main_row = main_row,
    r_group1 = r1, r_group2 = r2, delta_r = delta_r, route = route,
    stringsAsFactors = FALSE
  )

  if (i %% 50 == 0) cat(sprintf("  [%d/%d] processed\n", i, attempted))
}

# ---------------------------- Save outputs -----------------------------------
results_df <- dplyr::bind_rows(results)
names(results_df)[names(results_df) == "r_group1"] <- paste0("r_", GROUP_LEVELS[1])
names(results_df)[names(results_df) == "r_group2"] <- paste0("r_", GROUP_LEVELS[2])

write.csv(results_df,
          file.path(OUT_ROOT, sprintf("%s_main_smorf_coexpression_results.csv", PREFIX)),
          row.names = FALSE)

global_lines_df <- dplyr::bind_rows(global_lines)
saveRDS(global_lines_df, file.path(OUT_ROOT, sprintf("%s_global_lines.rds", PREFIX)))

cat("\nRun summary:\n")
cat(sprintf("  Attempted:           %d\n", attempted))
cat(sprintf("  Tested/written:      %d\n", nrow(results_df)))
cat(sprintf("  Skipped (no row):    %d\n", n_skip_missing))
cat(sprintf("  Skipped (invalid):   %d\n", n_skip_invalid))
cat(sprintf("  Coupled figures:     %d\n", n_coupled))
cat(sprintf("  Non-coupled figures: %d\n", n_non_coupled))
cat(sprintf("Outputs in: %s\n", OUT_ROOT))
