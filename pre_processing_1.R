#! usr/R

##This script does the preprocessing of methylation data  to normalize and filter the data including sample QC, probes QC, removing sex chromosomes, SNPs probes, cross-reactive
#probes and correcting for cell compostition. No batch effect have been noted to require
#correction (see pre-processing.Rmd). Finally, the script calculates methylation values.
##input= Raw data .idat files
##output = mSet object (.RData format), normalized  M and beta values (csv format)

##Requirements
 library("minfi")
 library("IlluminaHumanMethylationEPICmanifest")

##choose p-value threshold for the signal quality
P = 0.01

##Name of the metharray sheet
metharray_sheet = "BEA18P208_MC_sample_sheet_20190321-2.csv"

##choose normalization method
## can be either: "preprocessFunnorm" or "preprocessQuantile"
Norm_method = "preprocessQuantile"

##Assign significant names to the columns (samples) that you extract from the targets table (check head(targets) to determine which column encodes better for your samples)
Sample_Name <- "Sample_Name"

##Load data (for more detailed descriptions, see pre-processing.Rmd)
dataDirectory <-  file.path("data/EWAS_raw_data/")
list.files(dataDirectory, recursive = TRUE)

##Reading the sample sheet in csv format (facilitates .idat files import). 
##Csv file contains one line of description per sample, including sample description and file path. 
##The basenames that are created are specific to each sample/.idat. Here we set targets then we read intensities (rgSet)
targets <- read.metharray.sheet(dataDirectory, pattern=metharray_sheet) 
rgSet <- read.metharray.exp(targets=targets)  

##Assign significant names to the columns (samples) that you extract from the targets table (check head(targets) to determine which column encodes better for your samples)
targets$ID <- targets$Sample_Name
sampleNames(rgSet) <- targets$ID

##QC: detection p-values of the signal quality 
detP <- detectionP(rgSet)
keep <- colMeans(detP) < P
rgSet <- rgSet[,keep]

# remove poor quality samples from targets data
targets <- targets[keep,]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]

# normalize the data; this results in a GenomicRatioSet object (Touleimat and Tost 2012)
if (Norm_method == "preprocessQuantile") {
  mSetSq <- preprocessQuantile(rgSet)
} 
if (Norm_method == "preprocessFunnorm") {
  mSetSq <- preprocessFunnorm(rgSet)
} 

##Filtering: Put probes in mSetSq and detP objects in the same order
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

##Filter probes failing in 1 or more sample (0.01), change threshold if wanted
keep <- rowSums(detP < P) == ncol(mSetSq) 

 ##Check number of probes to be kept/removed (to be output?)
table(keep)  
mSetSqFlt <- mSetSq[keep,]

##Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

##Save all necessary R objects for later use
write_csv(as.data.frame(mVals), "results/mVals.csv")
write_csv(as.data.frame(bVals), "results/bVals.csv")
save(mSetSqFlt, file = "results/mSetSqFlt.RData")
save(rgSet, file = "results/rgSet.RData")
