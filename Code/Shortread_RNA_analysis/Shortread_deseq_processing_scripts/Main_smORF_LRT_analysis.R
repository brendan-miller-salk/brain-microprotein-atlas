# =============================================================================
# Main-ORF / smORF co-expression LRT analysis
#
# For every smORF / main-ORF pair in the MASTER annotation table, this script:
#   1. Pulls batch-corrected VST expression for both members from the DESeq2
#      object produced by `RNA_differential_expression.R`.
#   2. Computes the Pearson correlation in each diagnosis group (default:
#      Control vs AD) and tests whether they differ via Fisher Z (psych::paired.r).
#   3. Fits two nested linear models per pair:
#         additive:    smORF ~ main_ORF + group
#         interaction: smORF ~ main_ORF * group
#      and reports the LRT p-value, F-statistic, and beta/SE for the AD term in
#      each model. The additive beta isolates a constant AD effect on smORF
#      expression after accounting for main-ORF expression; the interaction
#      term captures whether the slope itself differs by group.
#   4. Writes one summary CSV with all per-pair statistics
#      (`<prefix>_main_smorf_LRT_results.csv`).
#
# Generalizable: change CONFIG paths to switch between cohorts (ROSMAP, MSBB,
# etc.). The script consumes the `<prefix>_dds.rds` and
# `<prefix>_vsd_corrected.rds` produced upstream.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(psych)
})

# ---------------------------- CONFIG -----------------------------------------
# Paths are relative to DESEQ_DIR (DESeq2 outputs from RNA_differential_expression.R),
# DATA_DIR (annotation inputs), and OUT_DIR (this script's outputs). Override via env
# vars `AD_PAPER_DESEQ_DIR`, `AD_PAPER_DATA_DIR`, `AD_PAPER_LRT_OUT_DIR`.
DESEQ_DIR <- Sys.getenv("AD_PAPER_DESEQ_DIR",   unset = "./outputs/deseq")
DATA_DIR  <- Sys.getenv("AD_PAPER_DATA_DIR",    unset = "./data")
OUT_DIR   <- Sys.getenv("AD_PAPER_LRT_OUT_DIR", unset = "./outputs/main_smorf_LRT")

DATASET <- "ROSMAP"
PREFIX  <- tolower(DATASET)

DDS_PATH <- file.path(DESEQ_DIR, sprintf("%s_dds.rds", PREFIX))
VSD_PATH <- file.path(DESEQ_DIR, sprintf("%s_vsd_corrected.rds", PREFIX))

# MASTER annotation: must contain a smORF coordinate / symbol column and
# `RNA_Ensembl_Parent` giving the main-ORF Ensembl gene ID.
MASTER_PATH <- file.path(DATA_DIR, "Brain_Microproteins_Discovery_MASTER_dataset.csv")

GROUP_COL     <- "new_dx"
GROUP_LEVELS  <- c("Control", "AD")          # exactly 2 groups for Fisher Z + LRT
EXCLUDE_TYPES <- c("lncRNA", "psORF", "eORF")
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- Load inputs ------------------------------------
master <- read.csv(MASTER_PATH, stringsAsFactors = FALSE)
dds      <- readRDS(DDS_PATH)
vst_expr <- readRDS(VSD_PATH)
meta_df  <- as.data.frame(colData(dds))

stopifnot(GROUP_COL %in% colnames(meta_df))
stopifnot(length(GROUP_LEVELS) == 2)

g1_samples <- rownames(meta_df[meta_df[[GROUP_COL]] == GROUP_LEVELS[1], ])
g2_samples <- rownames(meta_df[meta_df[[GROUP_COL]] == GROUP_LEVELS[2], ])
n1 <- length(g1_samples); n2 <- length(g2_samples)
cat(sprintf("%s: %d, %s: %d\n", GROUP_LEVELS[1], n1, GROUP_LEVELS[2], n2))

# ---------------------------- Build smORF / main-ORF pairs -------------------
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
    smorf_id    = smorf_match,
    main_id     = main_id_raw,
    main_gene   = gene_name
  ) %>%
  dplyr::distinct()

cat(sprintf("Pairs to test: %d\n", nrow(pairs_df)))

# ---------------------------- Helpers ----------------------------------------
find_row <- function(id, expr_rownames) {
  if (id %in% expr_rownames) return(id)
  regex <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", id)
  regex <- gsub(":", ".", regex)
  m <- grep(regex, expr_rownames, value = TRUE)
  if (length(m) > 0) return(m[1])
  NA_character_
}

fit_lrt <- function(s, p, grp) {
  df <- data.frame(s = s, p = p, grp = factor(grp, levels = GROUP_LEVELS))
  m_null <- lm(s ~ p,           data = df)
  m_add  <- lm(s ~ p + grp,     data = df)
  m_int  <- lm(s ~ p * grp,     data = df)
  add_an <- anova(m_null, m_add)
  int_an <- anova(m_add,  m_int)
  add_co <- summary(m_add)$coefficients
  int_co <- summary(m_int)$coefficients
  add_term <- grep(paste0("^grp", GROUP_LEVELS[2], "$"), rownames(add_co), value = TRUE)
  int_term <- grep(paste0("^p:grp", GROUP_LEVELS[2], "$"), rownames(int_co), value = TRUE)
  list(
    add_F  = add_an$F[2],            add_p  = add_an$`Pr(>F)`[2],
    int_F  = int_an$F[2],            int_p  = int_an$`Pr(>F)`[2],
    add_b  = if (length(add_term)) add_co[add_term, "Estimate"]   else NA,
    add_se = if (length(add_term)) add_co[add_term, "Std. Error"] else NA,
    int_b  = if (length(int_term)) int_co[int_term, "Estimate"]   else NA,
    int_se = if (length(int_term)) int_co[int_term, "Std. Error"] else NA
  )
}

# ---------------------------- Main loop --------------------------------------
rn <- rownames(vst_expr)
results <- vector("list", nrow(pairs_df))

for (i in seq_len(nrow(pairs_df))) {
  smorf_raw <- pairs_df$smorf_id[i]
  main_raw  <- pairs_df$main_id[i]
  main_gene <- pairs_df$main_gene[i]

  smorf_row <- find_row(smorf_raw, rn)
  main_row  <- find_row(main_raw,  rn)
  if (is.na(smorf_row) || is.na(main_row)) next

  s1 <- as.numeric(vst_expr[smorf_row, g1_samples])
  p1 <- as.numeric(vst_expr[main_row,  g1_samples])
  s2 <- as.numeric(vst_expr[smorf_row, g2_samples])
  p2 <- as.numeric(vst_expr[main_row,  g2_samples])

  if (sum(stats::complete.cases(s1, p1)) < 3 ||
      sum(stats::complete.cases(s2, p2)) < 3) next

  r1 <- suppressWarnings(stats::cor(s1, p1, use = "complete.obs"))
  r2 <- suppressWarnings(stats::cor(s2, p2, use = "complete.obs"))
  if (is.na(r1) || is.na(r2)) next

  fz <- tryCatch(psych::paired.r(r1, r2, n = n1, n2 = n2, twotailed = TRUE),
                 error = function(e) list(z = NA, p = NA))
  lrt <- tryCatch(
    fit_lrt(c(s1, s2), c(p1, p2), c(rep(GROUP_LEVELS[1], n1), rep(GROUP_LEVELS[2], n2))),
    error = function(e) list(add_F = NA, add_p = NA, int_F = NA, int_p = NA,
                             add_b = NA, add_se = NA, int_b = NA, int_se = NA)
  )

  results[[i]] <- data.frame(
    smorf_id              = smorf_raw,
    main_id               = main_raw,
    main_gene_name        = main_gene,
    smorf_row             = smorf_row,
    main_row              = main_row,
    r_group1              = r1,
    r_group2              = r2,
    delta_r               = r2 - r1,
    n_group1              = n1,
    n_group2              = n2,
    fisher_z              = fz$z,
    fisher_p              = fz$p,
    lrt_additive_F        = lrt$add_F,
    lrt_additive_p        = lrt$add_p,
    lrt_interaction_F     = lrt$int_F,
    lrt_interaction_p     = lrt$int_p,
    additive_beta         = lrt$add_b,
    additive_beta_se      = lrt$add_se,
    interaction_beta      = lrt$int_b,
    interaction_beta_se   = lrt$int_se,
    stringsAsFactors = FALSE
  )

  if (i %% 50 == 0)
    cat(sprintf("  [%d/%d] processed\n", i, nrow(pairs_df)))
}

results_df <- dplyr::bind_rows(results) %>%
  dplyr::arrange(fisher_p)

# Group-name aware column rename for clarity in the saved table.
names(results_df)[names(results_df) == "r_group1"] <- paste0("r_", GROUP_LEVELS[1])
names(results_df)[names(results_df) == "r_group2"] <- paste0("r_", GROUP_LEVELS[2])
names(results_df)[names(results_df) == "n_group1"] <- paste0("n_", GROUP_LEVELS[1])
names(results_df)[names(results_df) == "n_group2"] <- paste0("n_", GROUP_LEVELS[2])

out_csv <- file.path(OUT_DIR, sprintf("%s_main_smorf_LRT_results.csv", PREFIX))
write.csv(results_df, out_csv, row.names = FALSE)
cat(sprintf("Wrote %d rows -> %s\n", nrow(results_df), out_csv))
