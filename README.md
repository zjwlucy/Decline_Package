# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * `___IID___`:                    individual ID
  * `___FID___`:                    family ID     (For unralted individuals, create a columns of the same number as "FID", e.g. a column of "1")
  * `___pre_fev1___`:               FEV1           
  * `___SNPs___`:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
  * `___timefactor_spiro___`:       time since baseline exam
  * `___age___`:                    time-varying age
  * `___smoking_status___`:         Time-varying smoking status - never, former, current  (also used as grouping variable for glmmkin)

### Baseline 
  * `___age_baseline___`:           baseline age
  * `___ht_baseline___`:            baseline height (in cm)
  * `___smoking_packyears_base___`: pack-years at baseline
  * `___sex___`:                    biological sex
 
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
