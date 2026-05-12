library(parallel)
library(doParallel)

# NORMALIZE RAW DATA
dat <- read.csv('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/round1/shortstop_proteogenomics_appended_raw_no_correction_intensity.csv')
# Set $X to rownames
rownames(dat) <- dat$X
# Remove $X
dat <- dat[, -1]

# Assuming merged_data is your merged dataframe
rename_columns <- function(colnames) {
  sapply(colnames, function(name) {
    gsub("^(b)(\\d)(\\.\\d+[NC]?)$", "\\10\\2\\3", name)
  })
}

# Apply the renaming function to the column names
colnames(dat) <- rename_columns(colnames(dat))

# Check the updated column names
colnames(dat)

datadir <- "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/shortstop_proteogenomics_appended_results/round1/"
numericMeta<-read.csv(file=paste0("/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/rosmap_metadata/rosmap_50batch_specimen_metadata_for_batch_correction.csv"),row.names=1,header=TRUE)
outputfigs <- outputtabs <- paste0(datadir)

## Perform Batch Correction using TAMPOR (central tendency two-way table median polish algorithm)
#    for algorithm description see https://www.github.com/edammer/TAMPOR
#    and ROSMAP DLPFC Proteomics Release method description (extended wiki) PDF
source(paste0('/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/tmt/R_scripts/',"TAMPOR_custom.R"))
TAMPORlist.GIShybrid<-TAMPOR(dat,numericMeta,noGIS=FALSE,useAllNonGIS=TRUE,GISchannels=c("126","131"),batchPrefixInSampleNames=TRUE,parallelThreads=8,outputSuffix="GIShybrid")
#Note missingness at <50% is always enforced, losing 2855 protein isoforms, in order to get recommended data cleaned of batch effects.

# Overwrite visualization with all 3 pages of recorded/captured output (RStudio graphics output panel must have been scaled sufficiently for capture):
pdf(paste0(outputfigs,"TAMPOR_smORFs_unique_PSMs.pdf"),width=18,height=12)
TAMPORlist.GIShybrid$MDSplots
TAMPORlist.GIShybrid$meanSDplots
TAMPORlist.GIShybrid$convergencePlots
dev.off()

## Remove GIS samples
numericMeta<-TAMPORlist.GIShybrid$traits
numericMeta<-numericMeta[-which(numericMeta$Emory.Coded.Dx.2017=="GIS"),]
cleanDat<-TAMPORlist.GIShybrid$cleanDat
cleanRelAbun<-TAMPORlist.GIShybrid$cleanRelAbun
cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))]
cleanRelAbun<-cleanRelAbun[,match(rownames(numericMeta),colnames(cleanRelAbun))]


write.csv(cleanRelAbun,paste0(outputtabs,"C1.median_polish_corrected_relative_reporter_abundance_SHORTSTOP_AND_PROTEO_ROUND1.csv"))
write.csv(cleanDat,paste0(outputtabs,"C2.median_polish_corrected_log2abundanceRatioCenteredOnMedianOfBatchMediansPerProtein_SHORTSTOP_AND_PROTEO_ROUND1.csv"))



