library(dplyr)

# Define the function to process .tsv files
process_and_merge_psm_counts <- function(files, output_csv) {
  # Initialize an empty list to store data frames
  df_list <- list()
  
  # Read each file and store the data frame in the list
  for (file in files) {
    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df_list[[file]] <- df
  }
  
  # Merge all data frames by "Index", keeping all rows
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Index", all = TRUE), df_list)
  
  # Create rownames as Index
  rownames(merged_df) <- merged_df$Index
  
  # Create a new column called NumberPSM_all that sums the NumberPSM column
  merged_df$NumberPSM_all <- rowSums(merged_df[, grepl("NumberPSM", names(merged_df))], na.rm = TRUE)
  
  # Create a new dataframe with Index and NumberPSM
  psm_counts <- merged_df[, c("Index", "NumberPSM_all")] %>%
    rename("gene_id" = "Index", "psm_count" = "NumberPSM_all")
  
  # Save psm_counts to a .csv file
  write.csv(psm_counts, output_csv, row.names = FALSE)
  
  # Remove columns that start with specific patterns
  columns_to_remove <- grep("^(Index|NumberPSM|Gene|MaxPepProb|ReferenceIntensity)", names(merged_df))
  processed_df <- merged_df[, -columns_to_remove]
  
  return(processed_df)
}

# Base directory
base_dir <- "/scratch1/brendajm/tmt_rosmap/round2"
output_csv <- file.path(base_dir, "shortstop_proteogenomics_appended_merged_psm_counts.csv")

# List to store all file paths
files <- list()

# Loop through each subdirectory to collect file paths
for (i in 1:50) {
  subdir <- file.path(base_dir, paste0("b", i), "shortstop_proteogenomics_appended_results_cpm05/tmt-report")
  file1 <- file.path(subdir, "abundance_protein_MD.tsv")
  
  if (file.exists(file1)) {
    files <- c(files, file1)
  } else {
    message("File not found: ", file1)
  }
}

# Process and merge the files
merged_data <- process_and_merge_psm_counts(files, output_csv)

# Print the first few rows of the merged dataframe
head(merged_data)

write.csv(merged_data, file = "/scratch1/brendajm/tmt_rosmap/round2/shortstop_proteogenomics_rescored_appended_results/shortstop_proteogenomics_appended_raw_no_correction_intensity.csv", row.names = TRUE)

#Rscript process_batch_reports_shortstop_proteogenomics_appended.R
