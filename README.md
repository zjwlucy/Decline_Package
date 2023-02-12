# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * __`IID`__:                    individual ID
  * __`FID`__:                    family ID     (For unralted individuals, create a columns of the same number as "FID", e.g. a column of "1")
  * __`pre_fev1`__:               FEV1           
  * __`SNPs`__:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
  * __`timefactor_spiro`__:       time since baseline exam
  * __`age`__:                    time-varying age
  * __`smoking_status`__:         Time-varying smoking status - never, former, current  (also used as grouping variable for glmmkin)

### Baseline 
  * __`age_baseline`___:           baseline age
  * ___`ht_baseline`___:            baseline height (in cm)
  * ___`smoking_packyears_base`___: pack-years at baseline
  * ___`sex`___:                    biological sex
 
### Other cohort-specific variables
  * race
  * PCs (PC1, PC2, ... ...) 
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
