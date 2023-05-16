# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * __`IID`__:                    individual ID (numeric)
  * __`FID`__:                    family ID (numeric); For unrelated individuals, create a column of the same number as "FID", e.g. a column of 1.
  * __`pre_fev1`__:               FEV1           
  * __`SNPs`__:                   SNP information, column name MUST starts with the lower case "rs", e.g. rs507211, rsChrPosRefAlt
  * __`timefactor_spiro`__:       time since the baseline exam (in YEARS)
  * __`age`__:                    time-varying age
  * __`smoking_status`__:         time-varying smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin

### Baseline 
  * __`age_baseline`__:           baseline age
  * __`ht_baseline`__:            baseline height (in cm)
  * __`smoking_packyears_base`__: pack-years at baseline
  * __`sex`__:                    biological sex (female=0, male=1)
 
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
  * ggplot2





## For analysis, use "Analysis_MAIN.R"
