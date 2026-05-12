library(dplyr)
library(limma)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(reshape2)
library(doParallel)
library(parallel)
library(ggrepel)
#MAKE SURE YOU ADJUST DIMENSIONS OF CORRECTION

tmt_meta_file = "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/rosmap_metadata/combine_rounds_metadata.csv"
clinical_file = "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/rosmap_metadata/custom_dx_for_analyses.csv"
tmt_quant_file = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/C2.median_polish_corrected_log2abundanceRatioCenteredOnMedianOfBatchMediansPerProtein_SHORTSTOP_AND_PROTEO_COMBINEDROUNDS.csv'

process_tmt_data <- function(tmt_meta_file, clinical_file, tmt_quant_file) {
  # Load and rename columns in TMT metadata
  tmt_meta <- read.csv(tmt_meta_file)

  # Samples to remove from pioneering round of regression
  samples_to_remove <- c(
    "round1b22.129N", 
    "round1b27.128C", 
    "round1b36.127N", 
    "round1b41.127N", 
    "round1b43.127C", 
    "round2b02.129C", 
    "round2b06.131C", 
    "round2b09.128N"
  )  
   exclude <- tmt_meta %>%
     select(batch.channel) %>%
     as.vector()
   
   exclude$batch.channel <- unique(c(exclude$batch.channel, samples_to_remove))
  
  # Load and process the clinical data
  clinical <- read.csv(clinical_file) %>%
    mutate(age_death = ifelse(age_death == "90+", 91, age_death)) %>%
    mutate(age_death = as.numeric(age_death))
  
  
  # Load and process the TMT quantitation data
  tmt_quant <- read.csv(tmt_quant_file) %>%
    dplyr::rename(gene_id = X) %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>%
    tidyr::pivot_wider(names_from = gene_id, values_from = value) %>%
    dplyr::rename(batch.channel = sample) %>%
    dplyr::mutate(batch.channel = gsub("^b0", "b", batch.channel)) %>% 
    filter(!batch.channel %in% exclude$batch.channel)
  
  tmt_covariables <- clinical %>%
    select(c(batch.channel, age_death, msex, pmi, new_dx, Batch)) %>%
    filter(!is.na(new_dx)) %>%
    mutate(batch.channel = paste0(Batch, batch.channel)) %>%
    select(-c("Batch"))
  
  #  Merge the TMT data with the metadata
  tmt_data <- merge(tmt_quant, tmt_covariables, by = "batch.channel") %>%
    mutate(
      pmi_missing = ifelse(is.na(pmi), 1, 0),
      pmi = ifelse(is.na(pmi), max(pmi, na.rm = TRUE), pmi)
    )
  
  tmt_data_t <- t(tmt_data) %>% as.data.frame()
  tmt_data_t$protein <- rownames(tmt_data_t)
  
  # Convert EmoryStrictDx.2019 to a factor (categorical variable)
  tmt_data$new_dx <- as.factor(tmt_data$new_dx)
  
  # Enumerate the Diagnosis values (assign numeric codes)
  tmt_data$new_dx_num <- as.numeric(tmt_data$new_dx)
  
  # Number of bootstrap iterations
  n_bootstrap <- 5
  
  # Create an empty list to store bootstrap results
  bootstrap_results <- list()
  
  # Iterate over each protein (columns 3 to 50)
  #i in 2:7444 (MAKE SURE YOU ADJUST)
  for (i in 2:7444) {
    # Create an empty matrix to store coefficients for each bootstrap iteration
    bootstrap_coefficients <- matrix(NA, nrow = n_bootstrap, ncol = length(coef(lm(tmt_data[, i] ~ age_death + msex + pmi + new_dx, data = tmt_data))))
    colnames(bootstrap_coefficients) <- names(coef(lm(tmt_data[, i] ~ age_death + msex + pmi + new_dx, data = tmt_data)))
    
    for (b in 1:n_bootstrap) {
      # Sample with replacement
      sample_indices <- sample(1:nrow(tmt_data), replace = TRUE)
      bootstrap_sample <- tmt_data[sample_indices, ]
      
      # Fit the linear model to the bootstrap sample
      model <- lm(bootstrap_sample[, i] ~ age_death + msex + pmi + new_dx, data = bootstrap_sample)
      
      # Store the coefficients
      bootstrap_coefficients[b, ] <- coef(model)
    }
    
    # Convert the matrix to a dataframe and add to the list
    bootstrap_results[[colnames(tmt_data)[i]]] <- as.data.frame(bootstrap_coefficients)
  }
  
  # Calculate the median coefficients for each protein
  median_coefficients <- lapply(bootstrap_results, function(df) apply(df, 2, median))
  
  # Convert the list of median coefficients to a dataframe
  median_coefficients_df <- bind_rows(median_coefficients, .id = "Protein")
  
  # Correct the initial tmt_data using the median coefficients
  corrected_tmt_data <- tmt_data
  
  # Iterate over each protein (columns 3 to 50)
  for (i in 2:7444) {
    # Extract the median coefficients for the current protein
    protein_name <- colnames(tmt_data)[i]
    median_coeffs <- median_coefficients_df %>% filter(Protein == protein_name)
    
    # Correct the protein values in the initial tmt_data
    corrected_tmt_data[, i] <- tmt_data[, i] - (median_coeffs$age_death * tmt_data$age_death +
                                                  median_coeffs$msex * tmt_data$msex +
                                                  median_coeffs$pmi * tmt_data$pmi)
  }
  
  # View the first few rows of the corrected tmt_data
  head(corrected_tmt_data)
  return(corrected_tmt_data)
  
}


corrected_tmt_data <- process_tmt_data(
  tmt_meta_file,
  clinical_file,
  tmt_quant_file)

results_fit <- corrected_tmt_data %>% t() 
results_fit <- as.matrix(results_fit)
colnames(results_fit) <- corrected_tmt_data$batch.channel
results_fit <- results_fit[-1,]  
dim(results_fit)
results_fit <- results_fit[-c(7444:7450),]
write.csv(results_fit, "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/regressed_out_covariates_TMT_microproteins_shortstop_and_proteo_combined_rounds_PCAcorrected.csv")

# Transpose corrected data and add rownames as a new column
corrected_tmt_data_t <- corrected_tmt_data %>%
  t() %>%
  as.data.frame() %>%
  mutate(protein = rownames(.))

# Set rownames of corrected_tmt_data and prepare data for PCA
rownames(corrected_tmt_data) <- corrected_tmt_data$batch.channel
tmt_for_pca <- corrected_tmt_data[, 1:7444] %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ as.numeric(.)))

# Count missing and infinite values
cat("Missing values:", sum(is.na(tmt_for_pca)), "\n")
cat("Infinite values:", sum(is.infinite(tmt_for_pca)), "\n")

# Replace missing and infinite values with column means
tmt_for_pca_cleaned <- tmt_for_pca %>%
  mutate(across(everything(), ~ ifelse(is.na(.) | is.infinite(.), mean(., na.rm = TRUE), .)))

# Remove zero-variance columns
# Step 1: Identify zero-variance columns (handle NA)
zero_variance_cols <- apply(tmt_for_pca_cleaned, 2, function(x) var(x, na.rm = TRUE) == 0)

# Step 2: Replace NA in zero_variance_cols with TRUE (treat as zero variance)
zero_variance_cols[is.na(zero_variance_cols)] <- TRUE

# Step 3: Remove zero-variance columns
tmt_for_pca_cleaned <- tmt_for_pca_cleaned[, !zero_variance_cols]

# Perform PCA
pca <- prcomp(tmt_for_pca_cleaned, scale. = TRUE)
cat("PCA Summary:\n")
print(summary(pca))

# Extract and process PCA scores
pca_scores <- pca$x
pca_z_scores <- scale(pca_scores)
pca_scores <- pca_scores[, 1:5]

# Identify extreme samples (z > 2.5 or z < -2.5) in PC1 or PC2
extreme_samples <- apply(pca_z_scores[, 1:2], 1, function(row) any(abs(row) > 3))
extreme_sample_indices <- which(extreme_samples)
extreme_sample_names <- rownames(tmt_for_pca_cleaned)[extreme_sample_indices]

# Print names of extreme samples
cat("Extreme samples:\n")
print(extreme_sample_names)

# Plot PCA with extreme samples highlighted
pca_df <- data.frame(PC1 = pca_scores[, 1], PC2 = pca_scores[, 2], 
                     sample = rownames(tmt_for_pca_cleaned),
                     extreme = extreme_samples)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = extreme)) +
  geom_text_repel(data = pca_df[extreme_samples,], aes(label = sample), box.padding = 0.5) +
  theme_minimal() +
  labs(title = "PCA Plot with Extreme Samples Highlighted", color = "Extreme")

# # Calculate correlations
# correlations <- cor(pc_meta_data)
# correlations_melted <- melt(correlations)
# 
# # Plot the correlations
# correlation_plot <- ggplot(correlations_melted, aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1, 1), space = "Lab", 
#                        name = "Correlation") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1)) +
#   coord_fixed()