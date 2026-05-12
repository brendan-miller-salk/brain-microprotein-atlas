library(dplyr)
library(parallel)
library(doParallel)


# NORMALIZE RAW DATA
dat1 <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/round1/C1.median_polish_corrected_relative_reporter_abundance_SHORTSTOP_AND_PROTEO_ROUND1.csv')
dim(dat1)
# Set $X to rownames
rownames(dat1) <- dat1$X
# Remove $X
dat1 <- dat1[, -1]
dim(dat1)

dat2 <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/round2/C1.median_polish_corrected_relative_reporter_abundance_SHORTSTOP_AND_PROTEO_ROUND2.csv')
rownames(dat2) <- dat2$X
dat2 <- dat2[, -1]
dim(dat2)


# Add "round1" prefix to dat1 column names
colnames(dat1) <- paste0("round1", colnames(dat1))

# Add "round2" prefix to dat2 column names
colnames(dat2) <- paste0("round2", colnames(dat2))

# Combine the data frames by row names (protein IDs)
library(tibble)
combined_df <- full_join(
  dat1 %>% rownames_to_column("Protein"),
  dat2 %>% rownames_to_column("Protein"),
  by = "Protein"
)

# Replace NA with 0 if needed
combined_df[is.na(combined_df)] <- NA

# Convert back to a data frame with Protein as rownames if necessary
rownames(combined_df) <- combined_df$Protein
combined_df$Protein <- NULL

# Check the combined data frame
head(combined_df)

# Assuming merged_data is your merged dataframe
rename_columns <- function(colnames) {
  sapply(colnames, function(name) {
    gsub("^(b)(\\d)(\\.\\d+[NC]?)$", "\\10\\2\\3", name)
  })
}

# Apply the renaming function to the column names
colnames(combined_df) <- rename_columns(colnames(combined_df))

# Check the updated column names
colnames(combined_df)

datadir <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_results/"
numericMeta<-read.csv(file=paste0("/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/rosmap_metadata/combine_rounds_metadata.csv"),header=TRUE) %>%
  filter(dx != "GIS")
# Apply gsub to every cell in the data frame
numericMeta$batch.channel <- paste0(numericMeta$Batch, numericMeta$batch.channel)
numericMeta$batch.channel
rownames(numericMeta) <- numericMeta$batch.channel



## Perform Batch Correction using TAMPOR (central tendency two-way table median polish algorithm)
#    for algorithm description see https://www.github.com/edammer/TAMPOR
#    and ROSMAP DLPFC Proteomics Release method description (extended wiki) PDF
source(paste0('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/normalization_scripts/',"TAMPOR_custom.R"))
TAMPORlist.GIShybrid<-TAMPOR(combined_df,numericMeta,noGIS=TRUE,useAllNonGIS=TRUE,GISchannels=c(""),batchPrefixInSampleNames=FALSE,parallelThreads=8,outputSuffix="GIShybrid")
#Note missingness at <50% is always enforced, losing 2855 protein isoforms, in order to get recommended data cleaned of batch effects.


datadir <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_results/"
outputfigs <- outputtabs <- paste0(datadir)

# Overwrite visualization with all 3 pages of recorded/captured output (RStudio graphics output panel must have been scaled sufficiently for capture):
pdf(paste0(outputfigs,"TAMPOR_smORFs_unique_PSMs.pdf"),width=18,height=12)
TAMPORlist.GIShybrid$MDSplots
TAMPORlist.GIShybrid$meanSDplots
TAMPORlist.GIShybrid$convergencePlots
dev.off()

## Remove GIS samples
numericMeta<-TAMPORlist.GIShybrid$traits
cleanDat<-TAMPORlist.GIShybrid$cleanDat
cleanRelAbun<-TAMPORlist.GIShybrid$cleanRelAbun
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]
cleanRelAbun<-cleanRelAbun[,match(rownames(numericMeta),colnames(cleanRelAbun))]

write.csv(cleanRelAbun,paste0(outputtabs,"C1.median_polish_corrected_relative_reporter_abundance_SHORTSTOP_COMBINEDROUNDS.csv"))
write.csv(cleanDat,paste0(outputtabs,"C2.median_polish_corrected_log2abundanceRatioCenteredOnMedianOfBatchMediansPerProtein_SHORTSTOP_COMBINEDROUNDS.csv"))



