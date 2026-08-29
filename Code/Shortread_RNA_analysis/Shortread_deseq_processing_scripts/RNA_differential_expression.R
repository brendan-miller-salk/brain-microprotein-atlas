# =============================================================================
# Bulk RNA-seq differential expression (DESeq2)
#
# Generalizable script: configure dataset-specific paths, metadata columns,
# and grouping variables in the CONFIG block below. Defaults reproduce the
# ROSMAP DLPFC analysis used in the manuscript. Set DATASET to "MSBB" (or
# point the path variables at MSBB inputs) to run the equivalent MSBB analysis.
#
# Outputs (saved to OUT_DIR):
#   <prefix>_dds.rds             DESeq2 object (full design)
#   <prefix>_vsd.rds              VST-transformed expression
#   <prefix>_vsd_corrected.rds    VST with batch + covariate effects removed
#                                 (used by downstream LRT / co-expression scripts)
#   <prefix>_deseq2_<A>_vs_<B>_results.csv   AD vs Control contrast table
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---------------------------- CONFIG -----------------------------------------
# All paths are relative to DATA_DIR / OUT_DIR so the script is portable. Set
# the env vars `AD_PAPER_DATA_DIR` / `AD_PAPER_OUT_DIR` to override defaults.
DATA_DIR <- Sys.getenv("AD_PAPER_DATA_DIR", unset = "./data")
OUT_DIR  <- Sys.getenv("AD_PAPER_OUT_DIR",  unset = "./outputs/deseq")

DATASET <- "ROSMAP"   # label used in output filenames; e.g., "ROSMAP" or "MSBB"
PREFIX  <- tolower(DATASET)

# Counts: a single matrix CSV (genes x samples). For ROSMAP we merge hg38
# (Ensembl) + hg19 (Salk smORF) count tables; provide the merged matrix here.
COUNTS_PATH <- file.path(DATA_DIR, sprintf("%s_combined_counts.csv", PREFIX))

# Sample metadata CSV. Must include SAMPLE_ID_COL plus the columns referenced
# in DESIGN_COVARIATES and GROUP_COL.
META_PATH     <- file.path(DATA_DIR, sprintf("%s_clinical_rnaseq_metadata.csv", PREFIX))
SAMPLE_ID_COL <- "specimenID"

# Differential expression design.
GROUP_COL          <- "new_dx"               # primary grouping variable
GROUP_LEVELS       <- c("Control", "AsymAD", "AD")
CONTRAST           <- c("AD", "Control")     # numerator, denominator
BATCH_COL          <- "sequencingBatch"      # used for limma::removeBatchEffect
DESIGN_COVARIATES  <- c("RIN", "sequencingBatch", "msex", "age_death", "pmi")
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- Load counts ------------------------------------
counts <- read.csv(COUNTS_PATH, row.names = 1, check.names = FALSE)
counts[is.na(counts)] <- 0
cat(sprintf("Counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# ---------------------------- Load + clean metadata --------------------------
meta <- read.csv(META_PATH, stringsAsFactors = FALSE)

shared <- intersect(colnames(counts), meta[[SAMPLE_ID_COL]])
counts <- counts[, shared]
meta   <- meta[match(shared, meta[[SAMPLE_ID_COL]]), , drop = FALSE]
rownames(meta) <- meta[[SAMPLE_ID_COL]]

# Coerce common covariate types and scale continuous ones to avoid collinearity.
for (col in intersect(c("age_death", "ageDeath"), colnames(meta))) {
  meta[[col]] <- ifelse(meta[[col]] == "90+", 91, meta[[col]])
  meta[[col]] <- as.numeric(meta[[col]])
}
for (col in intersect(c("RIN", "pmi"), colnames(meta))) {
  meta[[col]] <- as.numeric(meta[[col]])
}
for (col in intersect(c("msex", "sex", "sequencingBatch"), colnames(meta))) {
  meta[[col]] <- as.factor(meta[[col]])
}

required <- c(SAMPLE_ID_COL, GROUP_COL, DESIGN_COVARIATES)
required <- intersect(required, colnames(meta))
meta <- meta[stats::complete.cases(meta[, required, drop = FALSE]), , drop = FALSE]
for (col in intersect(c("RIN", "pmi", "age_death", "ageDeath"), DESIGN_COVARIATES)) {
  meta[[col]] <- as.numeric(scale(meta[[col]]))
}

# ---------------------------- Diagnosis classification -----------------------
# Cohort-specific rubric. ROSMAP uses cogdx + mmse, MSBB uses CDR.
classify_rosmap <- function(braak, cerad, mmse, cogdx) {
  cerad <- ifelse(cerad == 4, 0,
           ifelse(cerad == 3, 1,
           ifelse(cerad == 2, 2,
           ifelse(cerad == 1, 3, cerad))))
  if (any(is.na(c(cerad, braak, mmse, cogdx)))) return(NA)
  dementia <- (mmse < 24) | (cogdx >= 4)
  if (cerad <= 1 && braak <= 3 && !dementia) {
    if (braak == 3 && cerad != 0) return(NA)
    return("Control")
  } else if (cerad >= 1 && cerad <= 3 && braak >= 3 && !dementia) {
    return("AsymAD")
  } else if (cerad >= 2 && braak >= 3 && dementia) {
    return("AD")
  }
  NA
}

classify_msbb <- function(braak, cerad, cdr) {
  cerad <- ifelse(cerad == 4, 0,
           ifelse(cerad == 3, 1,
           ifelse(cerad == 2, 2,
           ifelse(cerad == 1, 3, cerad))))
  if (any(is.na(c(cerad, braak, cdr)))) return(NA)
  dementia <- cdr >= 1
  if (cerad >= 1 && cerad <= 3 && braak >= 3 && !dementia) return("AsymAD")
  if (cerad <= 3 && braak <= 3 && !dementia)              return("Control")
  if (cerad >= 1 && braak >= 3 && dementia)               return("AD")
  NA
}

# Only (re)classify if GROUP_COL is missing or all-NA — otherwise trust input.
if (!(GROUP_COL %in% colnames(meta)) || all(is.na(meta[[GROUP_COL]]))) {
  if (DATASET == "ROSMAP") {
    meta[[GROUP_COL]] <- mapply(classify_rosmap,
                                meta$braaksc, meta$ceradsc,
                                meta$cts_mmse30_lv, meta$cogdx)
  } else if (DATASET == "MSBB") {
    meta[[GROUP_COL]] <- mapply(classify_msbb, meta$Braak, meta$CERAD, meta$CDR)
  }
}
meta <- meta[!is.na(meta[[GROUP_COL]]), , drop = FALSE]
meta[[GROUP_COL]] <- factor(meta[[GROUP_COL]],
                            levels = intersect(GROUP_LEVELS, unique(meta[[GROUP_COL]])))

counts <- counts[, rownames(meta)]
cat(sprintf("Final: %d genes x %d samples\n", nrow(counts), ncol(counts)))
print(table(meta[[GROUP_COL]]))

# ---------------------------- DESeq2 -----------------------------------------
design_formula <- as.formula(paste("~", paste(c(DESIGN_COVARIATES, GROUP_COL), collapse = " + ")))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
dds <- DESeq(dds, parallel = TRUE)
saveRDS(dds, file.path(OUT_DIR, sprintf("%s_dds.rds", PREFIX)))

# ---------------------------- VST + batch-corrected VST ----------------------
# Saved for downstream LRT and main-ORF / smORF co-expression scripts.
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file.path(OUT_DIR, sprintf("%s_vsd.rds", PREFIX)))

covariate_terms <- setdiff(DESIGN_COVARIATES, BATCH_COL)
cov_formula <- as.formula(paste("~", paste(c(BATCH_COL, covariate_terms), collapse = " + ")))
cov_matrix  <- model.matrix(cov_formula, data = as.data.frame(colData(dds)))
vsd_corrected <- removeBatchEffect(
  assay(vsd),
  batch = colData(dds)[[BATCH_COL]],
  covariates = cov_matrix[, -1, drop = FALSE]
)
saveRDS(vsd_corrected, file.path(OUT_DIR, sprintf("%s_vsd_corrected.rds", PREFIX)))

# ---------------------------- Contrast results -------------------------------
res_df <- as.data.frame(results(dds, contrast = c(GROUP_COL, CONTRAST[1], CONTRAST[2])))
res_df$gene_id <- rownames(res_df)

res_df <- res_df %>%
  dplyr::select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, gene_id) %>%
  dplyr::rename_with(~ paste0(PREFIX, "RNA_", .), -gene_id)

# Map ENSG -> SYMBOL (leaves smORF / Salk IDs untouched).
ensg_clean <- sub("\\..*$", "", res_df$gene_id)
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensg_clean,
                                 column = "SYMBOL", keytype = "ENSEMBL",
                                 multiVals = "first")
res_df$gene_id <- ifelse(is.na(symbols), res_df$gene_id, symbols)

out_csv <- file.path(OUT_DIR, sprintf("%s_deseq2_%s_vs_%s_results.csv",
                                       PREFIX, CONTRAST[1], CONTRAST[2]))
write.csv(res_df, out_csv, row.names = FALSE)
cat(sprintf("Wrote %s\n", out_csv))
