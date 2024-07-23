###General script to correct for confounders: age and remove probes associated to age, BMI and
##input=Beta values.rdata
##output=Beta values.rdata 

library(dplyr)
library("readr")
library("minfi")
library(tibble)
library("ChAMP")
library(doParallel) 

registerDoParallel(makePSOCKcluster(3))

##Constant: p-value threshold 
t = 0.05

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
mVals <- read.csv("results/mVals_noXprob_noSNP.csv")
bVals <- read.csv("results/bVals_noXprob_noSNP.csv")
get(load("results/mSetSqFlt_noXprob_noSNP.RData"))
assays(mSetSqFlt, withDimnames = FALSE) [["M"]] <- mVals
assays(mSetSqFlt, withDimnames = FALSE) [["Beta"]] <- bVals

metadata <- read_csv("results/extensive_metadata.csv") %>%
            filter(sample_id %in% c(colnames(bVals))) 
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(mVals)), metadata$sample_id),]

###Set the metadata field for categorical confounder
CONFOUNDER <- metadata$age

###Removing probes with high correlation to a confounder
dmp.conf.sig <- dmpFinder(mVals, CONFOUNDER, type = "categorical") %>% 
                rownames_to_column("probe") %>%
                filter(qval < t) %>%
                dplyr::select(probe) 

###Remove probes related to sex, age and bmi
dmp.conf.sig <- as_tibble(rownames(mVals)) %>% 
  anti_join(dmp.conf.sig %>% 
              dplyr::select(probe),
            by = c("value" = "probe"))
mVals <- mVals[dmp.conf.sig$value,]
bVals <- bVals[dmp.conf.sig$value,]
mSetSqFlt <- mSetSqFlt[dmp.conf.sig$value,]
assays(mSetSqFlt, withDimnames = FALSE) [["M"]] <- mVals
assays(mSetSqFlt, withDimnames = FALSE) [["Beta"]] <- bVals

##Save the data
save(mSetSqFlt, file = "results/mSetSqFlt.filtered_probes.RData")
write_csv(bVals, "results/cbVals.filtered_probes.csv")
write_csv(mVals, "results/cmVals.filtered_probes.csv")

