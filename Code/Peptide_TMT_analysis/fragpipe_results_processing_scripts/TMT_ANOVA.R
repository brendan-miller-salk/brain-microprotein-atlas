# Load required libraries
library(dplyr)
library(tibble)
library(broom)
library(tidyr)
library(ggplot2)

# Read the normalized TMT intensity file
tmt <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/regressed_out_covariates_TMT_microproteins_shortstop_and_proteo_combined_rounds_PCAcorrected.csv')

# Read the diagnosis file and select relevant columns
dx <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/rosmap_metadata/custom_dx_for_analyses.csv') %>%
  mutate(batch.channel = paste0(Batch, batch.channel)) %>%
  dplyr::select(batch.channel, new_dx)

# Transform and merge data
transposed_df <- t(tmt) %>%
  as.data.frame() %>%
  setNames(tmt$X) %>%
  rownames_to_column("batch.channel") %>%
  dplyr::slice(-1) %>%
  mutate(across(-batch.channel, as.numeric)) %>%
  merge(dx, ., by = "batch.channel") %>%
  filter(!is.na(new_dx)) %>%
  mutate(diagnosis = new_dx) %>%
  dplyr::select(batch.channel, diagnosis, everything(), -new_dx)

Q5BLP8 <- transposed_df %>%
 dplyr::select(c("batch.channel", "diagnosis", "Q9HB66")) %>%
  dplyr::mutate(diagnosis = factor(diagnosis, levels = c("Control", "AsymAD", "AD")))


ggplot(Q5BLP8, aes(x = diagnosis, y = `Q9HB66`, fill = diagnosis)) +
  geom_boxplot(size = 1, color = 'black') +
  scale_fill_manual(values = c("Control" = "grey", "AsymAD" = "#e26868", "AD" = "#b70000")) +
  geom_point(size = 3, alpha = 0.5) +
  theme_classic() +
  labs(title = paste("AD-Regulated Microprotein"), y = "Log2 Norm Abundance", x = "") +
  theme(
    axis.text.x =element_text(size = 22, color = 'black', hjust = 1, angle = 45),
    axis.text.y =element_text(size = 20, color = 'black'),
    
    legend.position = "none",
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 25)
  )

# Load the ggpmisc package
library(ggpmisc)

# Create the plot with regression stats
ggplot(Q5BLP8, aes(x = `ESPRESSO_chr21_796_170-chr21:25881673-25896501_F:0_P:0_M`, y = `P05067-11`)) +
  geom_smooth(method = "lm", se = TRUE, color = 'red') +  # Linear regression line with confidence interval
  geom_point(size = 3, alpha = 0.5) +  # Data points
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               label.x.npc = "left", label.y.npc = 0.95, size = 5, parse = TRUE) +  # Regression stats (equation and R-squared)
  theme_classic() +  # Clean theme
  labs(title = paste("Correlation between APP and Unnamed Microprotein"), y = "APP (Log2)", x = "Unnamed Microprotein (Log2)") +  # Labels
  theme(
    axis.text.x = element_text(size = 15, color = 'black'),  # X-axis text
    axis.text.y = element_text(size = 15, color = 'black'),  # Y-axis text
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 16),  # Axis title font size
    plot.title = element_text(size = 16)  # Plot title font size
  )


sem <- function(x) {
  sd(x) / sqrt(length(x))
}

P02743 %>%
  group_by(diagnosis) %>%
  summarise(
    mean = mean(P02743),
    sem = sem(P02743),
    lower_ci = mean(P02743) - 1.96 * sem(P02743),
    upper_ci = mean(P02743) + 1.96 * sem(P02743))

# Function to perform t-test for each protein and apply Tukey's HSD post-hoc test
perform_tests <- function(df) {
  # Ensure the diagnosis column is a factor
  #df <- transposed_df
  #df <-df[,c(1:100)]
  df$diagnosis <- as.factor(df$diagnosis)
  
  # Initialize a list to store t-test results
  ttest_results_list <- list()
  
  # Initialize a list to store ANOVA results
  anova_results_list <- list()
  
  # Initialize a list to store Tukey HSD results
  tukey_results_list <- list()
  
  # Iterate over each protein column (starting from the 3rd column) for t-test
  df_cntl_ad <- df %>% filter(diagnosis %in% c("Control", "AD"))
  for (i in 3:ncol(df)) {
    protein_name <- colnames(df_cntl_ad)[i]
    
    # Fit t-test
    t_test <- t.test(df_cntl_ad[[i]] ~ diagnosis, data = df_cntl_ad, var.equal = TRUE)
    
    # Extract the t-test p-value and add protein name
    ttest_results <- broom::tidy(t_test) %>%
      mutate(protein = protein_name)
    
    ttest_results_list[[protein_name]] <- ttest_results
  }
  
  # Iterate over each protein column (starting from the 3rd column) for ANOVA and Tukey's HSD
  for (i in 3:ncol(df)) {
    protein_name <- colnames(df)[i]
    
    # Fit the one-way ANOVA model
    anova_model <- aov(df[[i]] ~ diagnosis, data = df)
    
    # Extract the ANOVA p-value and add protein name
    anova_results <- broom::tidy(anova(anova_model)) %>%
      filter(term == "diagnosis") %>%
      mutate(protein = protein_name)
    
    # Perform Tukey's HSD post-hoc test
    tukey_test <- TukeyHSD(anova_model, "diagnosis")
    
    # Tidy the Tukey HSD results
    tukey_results <- as.data.frame(tukey_test$diagnosis) %>%
      rownames_to_column("comparison") %>%
      mutate(protein = protein_name)
    
    # Store the results in the respective lists
    anova_results_list[[protein_name]] <- anova_results
    tukey_results_list[[protein_name]] <- tukey_results
  }
  
  # Combine all t-test results into a single dataframe
  ttest_results_df <- bind_rows(ttest_results_list)
  
  # Combine all ANOVA results into a single dataframe
  anova_results_df <- bind_rows(anova_results_list)
  
  # Combine all Tukey HSD results into a single dataframe
  tukey_results_df <- bind_rows(tukey_results_list)
  
  # Calculate Q values (adjusted p-values) using the Benjamini-Hochberg method for t-test p-values
  ttest_results_df <- ttest_results_df %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    arrange(q.value, p.value)
  
  anova_results_df <- anova_results_df %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    arrange(q.value, p.value)
  
  tukey_results_df <- tukey_results_df %>%
    mutate(q.value = p.adjust(`p adj`, method = "BH")) %>%
    arrange(q.value, `p adj`)
  
  return(list(ttest_results = ttest_results_df, anova_results = anova_results_df, tukey_results = tukey_results_df, df_cntl_ad =df_cntl_ad))
}

# Perform the tests on the full dataset
#test_df <- transposed_df[,c(1:100)]
results <- perform_tests(transposed_df)

# Extract the test results
## Get Gene IDS based on UniProt
library(biomaRt)
## get some UniProt IDs
library(UniProt.ws)
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
ttest_results <- results$ttest_results
ttest_results$protein_simple <- sub("-.*$", "", ttest_results$protein)
up <- getBM(c("hgnc_symbol","uniprot_gn_id"), "uniprot_gn_id", ttest_results$protein_simple, mart)
colnames(up)[] <- c("gene_name", "protein_simple")
ttest_results <- merge(ttest_results, up, all.x = T) 
ttest_results$gene_name <- ifelse(is.na(ttest_results$gene_name), ttest_results$protein, ttest_results$gene_name)
uniprot_sequences <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/reference_genomes/human_genomes/uniprotkb_proteome_UP000005640_2024_12_09.csv') %>%
  rename(protein = gene_id)
ttest_results_uniprot <- merge(ttest_results, uniprot_sequences, by = "protein", all.x = TRUE)
# Create the new dataframe with selected columns and renamed headers
ttest_results_uniprot <- ttest_results_uniprot %>%
  dplyr::select(
    TMT_protein = protein,
    TMT_log2fc = estimate,
    TMT_pvalue = p.value,
    TMT_qvalue = q.value,
    TMT_gene_name = gene_name,
    sequence
  )
View(ttest_results_uniprot)
write.csv(ttest_results_uniprot, '/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/ttest_control_vs_AD_tmt_shortstop_and_proteo_results_cpm05.csv', row.names = F)


anova_results <- results$anova_results
anova_results$protein_simple <- sub("-.*$", "", anova_results$protein)
up <- getBM(c("hgnc_symbol","uniprot_gn_id"), "uniprot_gn_id", anova_results$protein_simple, mart)
colnames(up)[] <- c("symbol", "protein_simple")
anova_results <- merge(anova_results, up, all.x = T) 

tukey_results <- results$tukey_results
tukey_results$protein_simple <- sub("-.*$", "", tukey_results$protein)
up <- getBM(c("hgnc_symbol","uniprot_gn_id"), "uniprot_gn_id", tukey_results$protein_simple, mart)
colnames(up)[] <- c("symbol", "protein_simple")
tukey_results <- merge(tukey_results, up, all.x = T) 

write.csv(anova_results, '/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/amova_tmt_shortstop_and_proteo_results_cpm05.csv', row.names = F)
write.csv(tukey_results, '/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/tukey_tmt_shortstop__and_proteo_results_cpm05.csv', row.names = F)

ggplot(transposed_df, aes(diagnosis, `ESPRESSO_chr21_796_170-chr21:25881673-25896501_F:0_P:0_M`)) +
  geom_boxplot()



# Flatten the entire abundance matrix
all_abundance_values <- as.vector(as.matrix(transposed_df[ , 3:ncol(transposed_df)]))

# Convert all values to percentiles
percentile_values <- pnorm(all_abundance_values)

# Reshape back to the original structure
percentile_matrix <- matrix(percentile_values, 
                            nrow = nrow(transposed_df[ , 3:ncol(transposed_df)]), 
                            ncol = ncol(transposed_df[ , 3:ncol(transposed_df)]))

# Replace the original abundance values with the percentiles
transposed_df[ , 3:ncol(transposed_df)] <- percentile_matrix


write.csv(transposed_df, '/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/transposed_intensities_shortstop_and_proteo_percentiles_cpm05.csv')
table(transposed_df$diagnosis)
