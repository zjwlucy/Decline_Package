# code for fitting models

##############################################################  
## ------------------------------------------------------------ 
## data
## https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#
## covariates:

## IID:                    individual ID
## FID:                    family ID  
## pre_fev1:               FEV1           
## SNPs:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
## timefactor_spiro:       time since baseline exam
## age:                    time-varying age
## smoking_status:         Time-varying - never, former, current  (also used as grouping variable for glmmkin)

## Baseline 
 # age_baseline:           baseline age
 # ht_baseline:            baseline height (in cm)
 # smoking_packyears_base: pack-years at baseline
 # sex:                    biological sex
 
## other cohort-specific variables.
 # race
 # PCs 
 # equipchange ......


## kinship matrix (for related data):   both row names and column names MUST be IID      
      

## required package
 # geepack
 # GMMAT
 # dplyr
 # readxl
 # kinship2 (optional) for kinship matrix


  

##############################################################  
#                 0. Input Needed from users                 #
##############################################################

  d    <-   # read the dataset
  kmat <-   # kinship matrix (set to NULL if individuals are unrelated)
  
  covars_cohort  <-             # cohort specific covariates: e.g. PCs, equipmentchange......  (set to NULL if no additional covariates)
  data_multiRace <-             # multiple races? (TRUE/FALSE)
  data_related   <-             # For GEE, dataset with related individuals or not?  (TRUE/FALSE)     


############################################################## 








##############################################################  
#           1. preparation/checking (must run)               #
##############################################################
  source("f_data.R")
 
  covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace)
  covars_for_eq <- covars$covars_for_eq        # use for fitting models
  covars_common <- covars$covars_common        # use for checking data


# checking and cleaning the data
  d1 <- check_data(d, covars2=covars_common)   
  d2 <- d_slope(d1)                           # Full slope data
  d3 <- d_slope(d1,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
 




##############################################################  
#                  2. fit models separately                  #
##############################################################


## -----------------------------------------------------------
## example code: SNP rs507211


## LME
## fit ALL lme models for one SNP:
   fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat,          covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat)
   
## fit model 3 and model 6 (LME) for one SNP:   
   fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat[c(3,6),], covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat) 





## GEE
## fit ALL GEE models for one SNP
   fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee,           covars_additional=covars_for_eq, snpi="rs507211", related=data_related)

## fit model 1 and 3 for one SNP
   fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee[c(1,3), ], covars_additional=covars_for_eq, snpi="rs507211", related=data_related)






##############################################################  
#                       3. fit ALL models                    #
##############################################################
## -----------------------------------------------------------

## Note: fitting ALL models for ALL SNPs may take a long time
   source("f_fit.R")


## fit ALL models for ALL SNPs but do not save the output
   output <- f_fit(saveOutput=FALSE)

 

## fit ALL models for ALL SNPs and save the output to the excel file "summary_2023.xlsx"
   output <- f_fit(saveOutput=TRUE)
   






