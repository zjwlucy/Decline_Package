# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * IID:                    individual ID
  * FID:                    family ID     (For unralted individuals, create a columns of the same number as "FID", e.g. a column of "1")
  * pre_fev1:               FEV1           
  * SNPs:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
  * timefactor_spiro:       time since baseline exam
  * age:                    time-varying age
  * smoking_status:         Time-varying smoking status - never, former, current  (also used as grouping variable for glmmkin)

### Baseline 
  * age_baseline:           baseline age
  * ht_baseline:            baseline height (in cm)
  * smoking_packyears_base: pack-years at baseline
  * sex:                    biological sex
 
### Other cohort-specific variables.
  * race
  * PCs 
  * equipchange ......  


  * Kinship matrix (for related data):   both row names and column names MUST be IID      
      




## R version: R/4.0.2

### Required package
  * geepack
  * GMMAT
  * dplyr
  * readxl
  * kinship2 (optional) for kinship matrix






## For analysis, use "Analysis_MAIN.R"
