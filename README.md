# Decline_Package
Analysis plan on:  https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#



## Columns for data set 
  * __`IID`__:                    individual ID (numeric)
  * __`FID`__:                    family ID (numeric); For unrelated individuals, create a column of the same number as "FID", e.g. a column of 1
  * __`pre_fev1`__:               FEV1           
  * __`SNPs`__:                   SNP information, column name MUST starts with the lower case "rs", e.g. rs507211, rsChrPosRefAlt
  * __`timefactor_spiro`__:       time since the baseline exam (in YEARS). At baseline, timefactor_spiro=0
  * __`age`__:                    time-varying age
  * __`smoking_status`__:         time-varying smoking status (never=0, former=1, current=2); 

### Baseline 
  * __`age_baseline`__:           baseline age
  * __`ht_baseline`__:            baseline height (in cm)
  * __`smoking_packyears_base`__: pack-years at baseline
  * __`sex`__:                    biological sex (female=0, male=1)
  * __`smoking_status_base`__:    baseline smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin.
     
### Other cohort-specific variables
  * race (must be categorical)
  * PCs (PC1, PC2, ... ...) 
  * equipchange ......  

### Variables for summary not for analysis
  * __`pre_fev1fvc`__:   ratio of fev1 and fvc (fev1/fvc)
  * __`fev1_pp`__:       fev1 percent predicted


###  

  * Kinship matrix (for related data):   both row names and column names MUST be IID      
      



## Required packages
If current version does not work, you can try previous version.
Both versions should give the SAME output

### R version: R/4.2.1     (current version)
   * geepack_1.3.9
   * GMMAT_1.4.0 
   * dplyr_1.1.2  
   * readxl_1.4.3
   * ggplot2_3.4.2     
   * kinship2_1.9.6 (optional) for kinship matrix 

### R version: R/4.0.2  (previous version)
   * geepack_1.3-2
   * GMMAT_1.3.1
   * dplyr_1.0.2
   * readxl_1.3.1
   * ggplot2_3.3.2
   * kinship2 (optional) 

## For analysis, use "Analysis_MAIN.R"
