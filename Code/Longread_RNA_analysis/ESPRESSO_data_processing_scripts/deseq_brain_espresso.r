# === Libraries ===
library(DESeq2)
library(dplyr)
library(tibble)
library(plotly)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(viridisLite)

# === Set Working Directory ===
setwd('your/project/path/here') # <-- Replace with relative path or project root

# === Load Count Data ===
espresso_counts <- read.csv('data/counts.csv', row.names = 1, check.names = FALSE) %>%
  dplyr::select(-c(1:2)) %>%
  mutate_all(as.integer) %>%
  replace(is.na(.), 0)

# === Load Filtered CPM Data ===
espresso_cpm05 <- read.csv('data/brain_espresso_medianCPM05.csv', row.names = 1, check.names = FALSE)
rownames(espresso_cpm05) <- gsub(':', '_', rownames(espresso_cpm05))

# === Filter Counts to CPM Transcripts ===
espresso_counts <- espresso_counts[rownames(espresso_cpm05), ]

# === Load Sample Metadata ===
sampleTable <- read.csv('data/nanopore_metadata.csv') %>%
  filter(!is.na(sample_id)) %>%
  arrange(match(sample_id, colnames(espresso_cpm05))) %>%
  mutate(condition = factor(condition))

rownames(sampleTable) <- sampleTable$sample_id

# === Run DESeq2 ===
dds <- DESeqDataSetFromMatrix(countData = espresso_counts, colData = sampleTable, design = ~ condition + sex)
dds <- DESeq(dds, parallel = TRUE)

# === Extract Results ===
res <- results(dds, name = "condition_CT_vs_AD")
res$log2FoldChange <- -res$log2FoldChange
res$gene_id <- rownames(res)

sig_genes <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

# === PCA and Heatmap ===
vsd <- varianceStabilizingTransformation(dds)
mat <- assay(vsd)[head(order(res$padj), nrow(sig_genes)), ]
mat <- mat - rowMeans(mat)

df <- as.data.frame(colData(vsd)[, c("condition", "sex")])
colnames(df)[1] <- "condition"
df$condition <- factor(df$condition)

pheatmap(mat,
         annotation_col = df,
         col = colorRampPalette(turbo(20))(20),
         scale = 'row',
         show_colnames = FALSE,
         show_rownames = FALSE)

# === PCA Plot ===
plotPCA(vsd, intgroup = c("condition", "sex")) +
  geom_point(aes(shape = factor(condition)), size = 10, alpha = 0.8) +
  labs(title = "PCA", x = "PC1 (44%)", y = "PC2 (10%)") +
  theme(
    axis.text = element_text(size = 15, color = 'black'),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 30),
    plot.subtitle = element_text(size = 18, hjust = 0, face = "italic"),
    plot.caption = element_text(size = 18, hjust = 0.5, face = "italic"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.key.size = unit(1.5, 'lines')
  )

# === Volcano Plot ===
res$Significance <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significance), alpha = 0.8) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# === Plot Counts for a Specific Transcript ===
plotCounts(dds,
           gene = "ESPRESSO_chr10_4020_24",
           intgroup = "condition",
           returnData = FALSE)

# === Export DESeq2 Results ===
res$transcript_id <- rownames(res)
write.csv(res, 'output/deseq_results.csv', row.names = FALSE)