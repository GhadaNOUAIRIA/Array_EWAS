#! usr/R
##This script reads M-values and Beta values of methylation arrays and corrects for cell composition (samples from whole blood)
##Use after visualizing cell composition in the different samples and deciding if you want to correct for that effect
##Input = Beta values
##Output = M-values and beta values corrected for cell commpositions

libraray("readr")
library("ChAMP")
library(lumi)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
read_csv("results/bVals.csv")

##Correct for cell composition using ChAMP 
bVals <- (champ.refbase(beta=bVals, arraytype="EPIC"))$CorrectedBeta
mVals <- beta2m(bVals)

##Save all necessary R objects for later use
write_csv(mVals, "results/cmVals.cell_comp.csv")
write_csv(bVals, "results/cbVals.cell_comp.csv")