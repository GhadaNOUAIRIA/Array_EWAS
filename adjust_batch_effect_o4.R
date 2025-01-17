#! usr/R
##This script reads M-values and Beta values of methylation arrays and corrects for different natural confounders (batch effects can be included as well)
##Use after visualizing different confounders and deciding about the one that need to be corrected for
##input = beta values, metadata with samples column under the name "sample_id", and names of the columns from metadata that you want to adjust the data for.
##output = beta values and M values adjusted for the confounders

library("ChAMP")
library(lumi)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
bVals <- read_csv("results/cbVals.cell_comp.csv")
metadata <- read_csv("results/extensive_metadata.csv")

###Define the confounders you want to adjust your data for
confounders = c("bmi_cat", "age")
  
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(bVals)), metadata$sample_id),]

##Correct for bmi and age using ChAMP 
bVals <- champ.runCombat(beta = as.data.frame(bVals),
                         pd = as.data.frame(metadata),
                         variablename = "group",
                         batchname = confounders,
                         logitTrans = TRUE) ####Change to FALSE if you are using M-values
###Log transformation to obtain M-values
mVals <- beta2m(bVals)

##Save all necessary R objects for later use
write.csv(mVals, "results/mVals.cell_comp.cor_bmi_age.csv")
write.csv(bVals, "results/bVals.cell_comp.cor_bmi_age.csv")
