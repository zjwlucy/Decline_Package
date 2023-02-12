# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * __IID__:                    individual ID
  * __FID__:                    family ID     (For unralted individuals, create a columns of the same number as "FID", e.g. a column of "1")
  * __pre_fev1__:               FEV1           
  * __SNPs__:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
  * __timefactor_spiro__:       time since baseline exam
  * __age__:                    time-varying age
  * __smoking_status__:         Time-varying smoking status - never, former, current  (also used as grouping variable for glmmkin)

### Baseline 
  * __age_baseline__:           baseline age
  * __ht_baseline__:            baseline height (in cm)
  * __smoking_packyears_base__: pack-years at baseline
  * __sex__:                    biological sex
 
### Other cohort-specific variables
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
