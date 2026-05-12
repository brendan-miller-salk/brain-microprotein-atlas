# ==============================================================================
# scRNAseq_summary.R
# Purpose:
#   - Build cell-type–resolved heatmaps for significant gene counts, mean log2FC,
#     and log10(PSM) using ROSMAP snRNA-seq results joined to microprotein
#     (dark proteome) annotations.
#   - Build per–cell-type volcano plots (Unreviewed only) colored by log10(PSM).
#   - Export individual figures and a combined heatmap+volcano panel.
#   - Create tidy tables summarizing volcano data by (sub)cell type.
#
# Notes:
#   - Computation and filtering logic matches the provided script.
#   - Style and organization improved for "GitHub quality".
#   - Uses relative paths for portability across systems.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(readr)
  library(scales)
  library(patchwork)
  library(rlang)
  library(tidyverse)
  library(glue)
  library(ComplexUpset)
})

# ------------------------------------------------------------------------------
# Config / Paths
# ------------------------------------------------------------------------------
OUT_DIR <- "../../Results/scRNA_Enrichment"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Theme helpers (visual only; does not change logic)
# ------------------------------------------------------------------------------
theme_min_publication <- function() {
  theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title       = element_text(size = 15, hjust = 0.5),
      axis.text        = element_text(color = "black")
    )
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# Check if required data files exist
if (!file.exists("../data/combined_gpath_results.csv")) {
  stop("ERROR: combined_gpath_results.csv not found in ../data/\n",
       "Please download this file from:\n",
       "https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/tree/main/Results\n",
       "and place it in the Code/data/ directory.")
}

if (!file.exists("../data/ac_list_mapping.csv")) {
  stop("ERROR: ac_list_mapping.csv not found in ../data/\n",
       "Please download this file from:\n", 
       "https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/tree/main/Results\n",
       "and place it in the Code/data/ directory.")
}

gpath <- read_csv("../data/combined_gpath_results.csv")
ac_map <- read_csv("../data/ac_list_mapping.csv")

dark_proteome <- read.csv("../data/microprotein_master.csv", stringsAsFactors = FALSE) %>%
  # ── Step 1: Remove contaminants ──────────────────────────────────────────────
  filter(!grepl("SHEEP|Peptide|Horse|BOVIN", gene_symbol, ignore.case = TRUE)) %>%
  # ── Step 2: Keep only relevant databases ─────────────────────────────────────
  filter(str_detect(Database, "TrEMBL|Salk|Swiss-Prot")) %>%
  # ── Step 3: Remap Database labels ────────────────────────────────────────────
  mutate(
    smorf_type = ifelse(Database == "TrEMBL", "TrEMBL", smorf_type),
    Database   = case_when(
      Database == "TrEMBL"     ~ "Salk",
      Database == "Swiss-Prot" ~ "Swiss-Prot-MP",
      TRUE                     ~ Database
    )
  ) %>%
  # ── Step 4: Build evidence flags ─────────────────────────────────────────────
  mutate(
    has_RiboSAM = (RiboCode == TRUE) &
      (shortstop_label %in% c("SAM-Secreted", "SAM-Intracellular")) &
      (!smorf_type %in% c("oORF")) &
      (!grepl("iso", smorf_type, ignore.case = TRUE)),

    has_LooseRiboSAM = (RiboCode == TRUE) &
      (shortstop_label %in% c("SAM-Secreted", "SAM-Intracellular")) &
      (!smorf_type %in% c("Iso")),

    DDA_evidence = total_unique_spectral_counts > 0,

    DIA_evidence = has_LooseRiboSAM &
      !is.na(Global.PG.Q.Value) & (Global.PG.Q.Value <= 0.01) &
      (Proteotypic == 1),

    has_MS = DDA_evidence | DIA_evidence
  ) %>%
  # ── Step 5: Filter to microproteins ──────────────────────────────────────────
  filter(protein_class_length == "Microprotein") %>%
  # ── Step 6: Split by database, apply Salk-specific filtering, dedup ──────────
  {
    mp_swiss <- filter(., Database == "Swiss-Prot-MP") %>%
      filter(!duplicated(sequence))

    mp_salk <- filter(., Database == "Salk") %>%
      filter(has_RiboSAM | DDA_evidence | DIA_evidence) %>%
      arrange(desc(has_MS)) %>%
      filter(!duplicated(sequence))

    bind_rows(mp_swiss, mp_salk)
  } %>%
  # ── Remap Database to Unreviewed/Swiss-Prot for downstream joins ───────────
  mutate(
    Database = case_when(
      Database == "Salk"         ~ "Unreviewed",
      Database == "Swiss-Prot-MP" ~ "Swiss-Prot",
      TRUE                        ~ Database
    )
  )

# Quick sanity check
print(table(dark_proteome$protein_class_length, dark_proteome$Database))

# ------------------------------------------------------------------------------
# Map AC gene names -> Ensembl IDs; classify cell types
# ------------------------------------------------------------------------------
gpath_clean <- gpath %>%
  left_join(ac_map, by = c("gene" = "gene_name")) %>%
  mutate(
    gene = if_else(!is.na(gene_id), gene_id, gene),
    cell_type_general = case_when(
      grepl("^Exc", celltype)                                           ~ "Exc Neur",
      grepl("^Inh", celltype)                                           ~ "Inh Neur",
      grepl("^Ast", celltype)                                           ~ "Astr",
      celltype %in% c("Mic","MicMKI67","MicP2RY12","MicTPT1","CAMs","Tcells") ~ "Immune",
      celltype %in% c("Oli","OPC")                                      ~ "Oli",
      celltype %in% c("End","Per","SMC","Fib")                         ~ "Vas",
      TRUE                                                              ~ "Other"
    )
  ) %>%
  select(-gene_id)

# ------------------------------------------------------------------------------
# Join with dark proteome metadata
# ------------------------------------------------------------------------------
gpath_annotated <- gpath_clean %>%
  left_join(
    dark_proteome %>%
      select(
        gene_name, gene_id, sequence, smorf_type, Database,
        total_unique_spectral_counts, total_razor_spectral_counts
      ),
    by = c("gene" = "gene_name")
  )

# ------------------------------------------------------------------------------
# Prepare plot data (significant genes only)
# ------------------------------------------------------------------------------
plot_data <- gpath_annotated %>%
  filter(p_adj.glb < 0.05) %>%
  mutate(
    cell_type_general = trimws(cell_type_general),
    celltype          = trimws(celltype)
  )

# Write main summary output (drop genes with no dark-proteome match)
write.csv(plot_data %>% filter(!is.na(Database)) %>% select(-c("Unnamed: 0")), file.path(OUT_DIR, 'scRNA_Enrichment_summary.csv'), row.names = FALSE)

# ------------------------------------------------------------------------------
# Aggregate for heatmaps
# ------------------------------------------------------------------------------
heatmap_data <- plot_data %>%
  filter(p_adj.glb < 0.05) %>%
  group_by(cell_type_general, celltype, smorf_type) %>%
  summarise(
    sig_gene_count = n(),
    log2fc = mean(logFC, na.rm = TRUE),
    total_unique_spectral_counts = case_when(
      smorf_type == "SwissProt" | smorf_type == "TrEMBL" ~ mean(total_razor_spectral_counts,  na.rm = TRUE),
      TRUE                      ~ mean(total_unique_spectral_counts, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Place "SwissProt" first in smorf_type
heatmap_data$smorf_type <- factor(
  heatmap_data$smorf_type,
  levels = c("SwissProt", setdiff(unique(heatmap_data$smorf_type), "SwissProt"))
)

# ------------------------------------------------------------------------------
# Heatmap 1: # Significant Genes
# ------------------------------------------------------------------------------
heatmap_number_genes <- ggplot(
  heatmap_data %>% filter(!is.na(smorf_type), smorf_type != ""),
  aes(x = smorf_type, y = fct_rev(celltype), fill = sig_gene_count)
) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "gold", high = "#d73027",
    midpoint = 50, limits = c(0, 100), oob = scales::squish,
    name = "# Significant Genes"
  ) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  facet_grid(
    rows = vars(cell_type_general),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_min_publication() +
  theme(
    axis.text.x.top    = element_text(angle = 45, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom= element_blank(),
    axis.title.x.bottom= element_blank(),
    axis.text.y.right  = element_text(size = 7, hjust = 0),
    axis.text.y.left   = element_blank(),
    axis.ticks.y.left  = element_blank(),
    axis.ticks.y.right = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(size = 8, angle = 90, hjust = 0.5),
    panel.spacing.y    = unit(1, "lines"),
    legend.position    = "bottom"
  )

ggsave(
  filename = file.path(OUT_DIR, "cell_type_smorf_type_heatmap.pdf"),
  plot = heatmap_number_genes,
  width = 8, height = 10
)

# ------------------------------------------------------------------------------
# Heatmap 2: mean log2FC
# ------------------------------------------------------------------------------
heatmap_log2FC <- ggplot(
  heatmap_data %>% filter(!is.na(smorf_type), smorf_type != ""),
  aes(x = smorf_type, y = fct_rev(celltype), fill = log2fc)
) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish,
    name = "Log2FC"
  ) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  facet_grid(
    rows = vars(cell_type_general),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_min_publication() +
  theme(
    axis.text.x.top    = element_text(angle = 45, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom= element_blank(),
    axis.title.x.bottom= element_blank(),
    axis.text.y.right  = element_text(size = 7, hjust = 0),
    axis.text.y.left   = element_blank(),
    axis.ticks.y.left  = element_blank(),
    axis.ticks.y.right = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(size = 8, angle = 90, hjust = 0.5),
    panel.spacing.y    = unit(1, "lines"),
    legend.position    = "bottom"
  )

ggsave(
  filename = file.path(OUT_DIR, "heatmap_log2FC.pdf"),
  plot = heatmap_log2FC,
  width = 8, height = 10
)

# ------------------------------------------------------------------------------
# Heatmap 3: log10(PSM)
# ------------------------------------------------------------------------------
heatmap_PSM <- ggplot(
  heatmap_data %>% filter(!is.na(smorf_type), smorf_type != ""),
  aes(x = smorf_type, y = fct_rev(celltype), fill = log10(total_unique_spectral_counts+0.1))
) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0.5, limits = c(0, 2),
    oob = scales::squish,
    name = "Log10PSM"
  ) + 
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  facet_grid(
    rows = vars(cell_type_general),
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_min_publication() +
  theme(
    axis.text.x.top    = element_text(angle = 45, hjust = 0),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom= element_blank(),
    axis.title.x.bottom= element_blank(),
    axis.text.y.right  = element_text(size = 7, hjust = 0),
    axis.text.y.left   = element_blank(),
    axis.ticks.y.left  = element_blank(),
    axis.ticks.y.right = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(size = 8, angle = 90, hjust = 0.5),
    panel.spacing.y    = unit(1, "lines"),
    legend.position    = "bottom"
  )

ggsave(
  filename = file.path(OUT_DIR, "heatmap_PSM.pdf"),
  plot = heatmap_PSM,
  width = 8, height = 10
)

# ------------------------------------------------------------------------------
# Volcano plots per cell_type_general (Unreviewed only)
# ------------------------------------------------------------------------------
ordered_cell_types <- c("Astr", "Exc Neur", "Immune", "Inh Neur", "Oli", "Vas")

plot_list <- list()

for (ct in ordered_cell_types) {
  message("Plotting volcano for: ", ct)
  
  df_ct <- gpath_annotated %>%
    filter(cell_type_general == ct) %>%
    filter(!is.na(Database)) %>%
    arrange(factor(Database, levels = c("Unreviewed", "Swiss-Prot"))) %>%
    filter(Database == "Unreviewed")
  
  if (nrow(df_ct) == 0) next
  
  gg <- ggplot(
    df_ct %>% filter(!is.na(Database)),
    aes(x = logFC, y = -log10(p_adj.loc))
  ) +
    geom_point(
      aes(
        fill  = log10(total_unique_spectral_counts+0.01),
        color = log10(total_unique_spectral_counts+0.01),
        alpha = ifelse(p_adj.loc > 0.05, 0.1, 0.85)
      ),
      shape = 21, size = 3, stroke = 0.6
    ) +
    geom_hline(yintercept = -log10(0.05), color = "navy", linewidth = 1) +
    scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 1.5, limits = c(0, 3), oob = scales::squish,
      name = "Log10PSM"
    ) +
    scale_color_gradient2(
      low = "navy", mid = "grey", high = "dark red",
      midpoint = 1, limits = c(0, 3), name = "Log10PSM"
    ) +
    scale_alpha_identity() +
    guides(
      color = "none",
      alpha = "none",
      fill  = guide_legend(override.aes = list(size = 3, shape = 21, stroke = 0.6))
    ) +
    labs(
      title = ct,
      x = "Log2 Fold Change",
      y = "-Log10 Q-value",
      fill = NULL
    ) +
    scale_x_continuous(limits = c(-1.2, 1.2)) +
    theme_min_publication() +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 14),
      axis.title.x     = element_text(size = 12),
      axis.title.y     = element_text(size = 12),
      axis.text        = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
  
  plot_list[[ct]] <- gg
}

if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    filename = file.path(OUT_DIR, "volcano_all_celltypes.pdf"),
    plot = combined_plot,
    width = 12, height = 8
  )
}

# ------------------------------------------------------------------------------
# Build volcano data tables per cell_type_general
# ------------------------------------------------------------------------------
volcano_data_list <- list()
cell_types_general <- unique(gpath_annotated$cell_type_general)

for (ct in cell_types_general) {
  message("Processing data for: ", ct)
  
  df_ct <- gpath_annotated %>%
    filter(cell_type_general == ct) %>%
    filter(!is.na(Database)) %>%
    mutate(log10_padj = -log10(p_adj.loc)) %>%
    select(gene, logFC, p_adj.loc, log10_padj, Database, smorf_type,
           celltype, cell_type_general, sequence, total_unique_spectral_counts)
  
  volcano_data_list[[ct]] <- df_ct
}

# ------------------------------------------------------------------------------
# Wide tables by (sub)celltype: one row per sequence, columns = logFC/p_adj.loc
# ------------------------------------------------------------------------------
simplified_by_celltype <- list()
celltypes <- unique(gpath_annotated$celltype)

for (ct in celltypes) {
  message("Processing celltype: ", ct)
  
  df_ct <- gpath_annotated %>%
    filter(celltype == ct, !is.na(Database))
  
  if (nrow(df_ct) == 0) next
  
  ct_label <- unique(df_ct$celltype)[1]
  
  df_wide <- df_ct %>%
    mutate(
      !!sym(paste0(ct_label, "_logFC"))     := logFC,
      !!sym(paste0(ct_label, "_p_adj.loc")) := p_adj.loc
    ) %>%
    select(sequence, ends_with("_logFC"), ends_with("_p_adj.loc"))
  
  simplified_by_celltype[[ct]] <- df_wide
}

final_celltype_df <- bind_rows(simplified_by_celltype)

# ------------------------------------------------------------------------------
# Wide tables by general cell type
# ------------------------------------------------------------------------------
simplified_by_general <- list()
generals <- unique(gpath_annotated$cell_type_general)

for (ct in generals) {
  message("Processing general cell type: ", ct)
  
  df_ct <- gpath_annotated %>%
    filter(cell_type_general == ct, !is.na(Database))
  
  if (nrow(df_ct) == 0) next
  
  ct_label <- gsub(" ", "_", unique(df_ct$cell_type_general)[1])
  
  df_wide <- df_ct %>%
    mutate(
      !!sym(paste0(ct_label, "_logFC"))     := logFC,
      !!sym(paste0(ct_label, "_p_adj.loc")) := p_adj.loc
    ) %>%
    select(sequence, ends_with("_logFC"), ends_with("_p_adj.loc"))
  
  simplified_by_general[[ct]] <- df_wide
}

final_general_df <- bind_rows(simplified_by_general)

# ------------------------------------------------------------------------------
# Merge back to master and flag scRNA significance anywhere (p_adj.loc < 0.05)
# ------------------------------------------------------------------------------
#dark_proteome_master <- read.csv("../data/microprotein_master.csv")

# if (nrow(final_general_df) > 0) {
#   dark_proteome_master <- merge(dark_proteome_master, final_general_df, by = "sequence")
#   dark_proteome_master$scrna_sign <- apply(
#     dark_proteome_master[, grep("_p_adj\\.loc$", colnames(dark_proteome_master))],
#     1,
#     function(pvals) any(pvals < 0.05, na.rm = TRUE)
#   )
# }

# Save final master dataset with scRNA annotations
#out_csv <- file.path(OUT_DIR, "Brain_Microproteins_Discovery_MASTER_dataset_scRNA.csv")
#write.csv(dark_proteome_master, out_csv, row.names = FALSE)
#message("Saved: ", out_csv)

#message("scRNAseq summary analysis completed!")

