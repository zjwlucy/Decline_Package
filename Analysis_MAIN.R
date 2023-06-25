# code for fitting models

##############################################################  
## ------------------------------------------------------------ 
## data
## https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#
## covariates:

## IID:                    individual ID (numeric)
## FID:                    family ID  (numeric)
## pre_fev1:               FEV1           
## SNPs:                   SNP information, column name MUST starts with lower case "rs", e.g. rs507211
## timefactor_spiro:       time since the baseline exam (in YEARS). At baseline, timefactor_spiro=0.
## age:                    time-varying age
## smoking_status:         time-varying smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin.

## Baseline 
 # age_baseline:           baseline age
 # ht_baseline:            baseline height (in cm)
 # smoking_packyears_base: pack-years at baseline
 # sex:                    biological sex (female=0, male=1)
 
## other cohort-specific variables.
 # race
 # PCs 
 # equipchange ......


## variables for summary not for analysis
## FEV1/FVC:               ratio of fev1 and fvc
## FEV1 % predicted:       fev1 percent predicted



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
  kmat <-   # kinship matrix (set to NULL if individuals are unrelated, e.g. kmat <- NULL)
  
  covars_cohort  <-             # cohort specific covariates: e.g. PCs, equipmentchange......  (set to NULL if no additional covariates)
  data_multiRace <-             # multiple races? (TRUE/FALSE)
  data_related   <-             # For GEE, dataset with related individuals or not?  (TRUE/FALSE)     
  cohortname     <-             # cohort name

## example:
#   d    <- read.csv("FHS_data.csv")  
#   kmat <- read.csv("FHS_kinship.csv")
   
#   covars_cohort  <-  c("PC1", "PC2", "equipchange")
#   data_multiRace <-  FALSE          
#   data_related   <-  TRUE   
#   cohortname     <-  "FHS"     

############################################################## 








##############################################################  
#           1. preparation/checking (must run)               #
##############################################################
  source("f_data.R")
 
  covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace)
  covars_for_eq <- covars$covars_for_eq        # use for fitting models
  covars_common <- covars$covars_common        # use for checking data





##############################################################  
#                  2. fit models separately                  #
##############################################################


## -----------------------------------------------------------
## example code: SNP rs507211

## checking and cleaning the data
   d1 <- check_data(d, covars2=covars_common, snpi="rs507211")   
   d2 <- d_slope(d1)                           # Full slope data
   d3 <- d_slope(d1,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
 

## LME
## fit ALL lme models for one SNP:
   out_gmmat1 <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat,          covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat)
   
## fit model 3 and model 6 (LME) for one SNP:   
   out_gmmat2 <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat[c(3,6),], covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat) 

## show coefficients for ALL variables
   out_gmmat3 <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat,          covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat, all_results=TRUE)


## fit base models for model 1 without SNPs
   out_gmmat4 <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, eqlist=l_gmmat[1,],      covars_additional=covars_for_eq, snpi="smoking_status", kmat=kmat) 
   

   out_gmmat1
   out_gmmat2
   out_gmmat3
   out_gmmat4



## GEE
## fit ALL GEE models for one SNP
   out_gee1 <- fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee,           covars_additional=covars_for_eq, snpi="rs507211", related=data_related)

## fit model 1 and 3 (GEE) for one SNP
   out_gee2 <- fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee[c(1,3), ], covars_additional=covars_for_eq, snpi="rs507211", related=data_related)

## show coefficients for ALL variables
   out_gee3 <- fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee,           covars_additional=covars_for_eq, snpi="rs507211", related=data_related, all_results=TRUE)


## fit base models for model 1 without SNPs
   out_gee4 <- fit_all_gee(dat_full=d1, dat_slope=d2, eqlist=l_gee[1, ],      covars_additional=covars_for_eq, snpi="smoking_status", related=data_related)


   out_gee1
   out_gee2
   out_gee3
   out_gee4






##############################################################  
#                       3. fit ALL models                    #
##############################################################
## -----------------------------------------------------------
## Make sure to run Section 1 before running the code below

## Note: fitting ALL models for ALL SNPs may take a long time
   source("f_fit.R")


## fit ALL the models for ALL SNPs but do not save the output
   output <- f_fit(saveOutput=FALSE)



## Model using SNP
## fit ALL the models for ALL SNPs and save SNP-related output to the excel file "summary_partial_2023.xlsx"
   output <- f_fit(saveOutput=TRUE, cohort=cohortname)
   
## fit ALL the models and save all the results to the excel file "summary_allresults_2023.xlsx"
   output <- f_fit(saveOutput=TRUE, allresults=TRUE, cohort=cohortname)
  
   
   
## Base model  
## fit ALL the BASE models without SNPs and save smoking_status-related output to the excel file "summary_partial_base_2023.xlsx"
   output <- f_fit(BaseModel=TRUE, voi="smoking_status", saveOutput=TRUE, cohort=cohortname)
  
## fit ALL the BASE models without SNPs and save smoking_status-related output to the excel file "summary_allresults_base_2023.xlsx"
   output <- f_fit(BaseModel=TRUE, voi="smoking_status", saveOutput=TRUE, allresults=TRUE, cohort=cohortname)




