# Array_EWAS
This is a collection of R scripts that handle methylation data from Illumina arrays. They perform pre-processing of the data, quality checks, confounder check, and find DMPs (differentially methylated positions) and DMRs (differentially methylated regions).

Order of the scripts:
pre_processing_1.R

##The following 2 scripts can come in any order
xreactive_probes_find_remove_2.R
remove_snp_probes_3.R 

##The following scripts are optional and can be used in any order
remove_sex_chromosomes_o1.R
remove_confounding_probes_o2.R
adjust_cell_composition_o3.R
adjust_batch_effect_o4.R


##The following 2 scripts are mandatory and should be run in the following order
find_dmp_4.R
find_dmr_5.R

##The following script can be optional
find_blocks_6.R
