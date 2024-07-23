#! usr/R

##This script is removing SNPs probes, Finally, the script calculates methylation values.
### input = normalized  and filtered mSet object 
##output = mSet object (.RData format), normalized  M and beta values (csv format)

##Requirements
library("minfi")

##Constants
MAF = 0

##Load data (for more detailed descriptions, see pre-processing.Rmd)
get(load("results/mSetSqFlt_noXprob.RData"))


##Remove probes with known SNPs, default is to remove all, but you can define a threshold of highest frequency of the minor allele
snps <- getSnpInfo(mSetSqFlt) 
mSetSqFlt <- addSnpInfo(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps=c("SBE","CpG", "Probe"), maf = MAF) 

##Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)


##Save all necessary R objects for later use
write_csv(as.data.frame(mVals), "results/mVals_noXprob_noSNP.csv")
write_csv(as.data.frame(bVals), "results/bVals_noXprob_noSNP.csv")
save(mSetSqFlt, file = "results/mSetSqFlt_noXprob_noSNP.RData")
