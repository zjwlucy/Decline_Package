# code for fitting models

##############################################################  
## ------------------------------------------------------------ 
## data
## https://docs.google.com/document/d/1ecNAiYQG7C4lmPHUhkvw3B1MMH9wfZikHFWfLRQGKtw/edit#
## covariates:

## IID:                    individual ID (numeric)
## FID:                    family ID  (numeric); For unrelated individuals, create a column of the same number as "FID", e.g. a column of 1
## pre_fev1:               FEV1 (in mL)         
## SNPs:                   SNP information, column name MUST starts with the lower case "rs", e.g. rs507211, rsChrPosRefAlt
## timefactor_spiro:       time since the baseline exam (in YEARS). At baseline, timefactor_spiro=0.
## age:                    time-varying age (in YEARS)
## smoking_status:         time-varying smoking status (never=0, former=1, current=2); 

## Baseline 
 # age_baseline:           baseline age
 # ht_baseline:            baseline height (in cm)
 # smoking_packyears_base: pack-years at baseline
 # sex:                    biological sex (female=0, male=1)
 # smoking_status_base:    baseline smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin.
 
## Other cohort-specific variables.
 # race (must be categorical)
 # PCs 
 # equipchange ......


## Variables for summary not for analysis
 # pre_fev1fvc:   ratio of fev1 and fvc (fev1/fvc)
 # fev1_pp:       fev1 percent predicted


## kinship matrix (for related data):   both row names and column names MUST be IID      
      



## ---------------------------------------------------------------
## Required packages
## If current version does not work, you can try previous version
## Both versions should give the SAME output


## R version: R/4.2.1     (current version)
   # geepack_1.3.9
   # GMMAT_1.4.0 
   # dplyr_1.1.2  
   # readxl_1.4.3
   # ggplot2_3.4.2     
   # kinship2_1.9.6 (optional) for kinship matrix 


## R version: R/4.0.2  (previous version)
   # geepack_1.3-2
   # GMMAT_1.3.1
   # dplyr_1.0.2
   # readxl_1.3.1
   # ggplot2_3.3.2
   # kinship2 (optional) 



 

  
   



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
#           (1). Preparation/checking (must run)             #
##############################################################
  source("f_data.R")
 
  covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace)
  covars_for_eq <- covars$covars_for_eq        # use for fitting models
  covars_common <- covars$covars_common        # use for checking data




##############################################################  
#                     (2). Fit ALL models                    #
##############################################################
## -----------------------------------------------------------
## Make sure to run Section (1) before running the code below
## DO NOT run section (3) before this section

## Note: 
 # Fitting ALL models for ALL SNPs may take a long time
   source("f_fit.R")


## fit ALL the models for ALL SNPs but do not save the output
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq, 
                   data_related=data_related,
                   saveOutput=FALSE)



## Model using SNP
## fit ALL the models for ALL SNPs and save SNP-related output to the csv file "Cohort_(GEE/GMMAT)_partial_2023.csv"
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq, 
                   data_related=data_related, 
                   saveOutput=TRUE, cohort=cohortname)
   
## fit ALL the models and save all the results to the csv file "Cohort_(GEE/GMMAT)_allresults_2023.csv"
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq,
                   data_related=data_related,
                   saveOutput=TRUE, cohort=cohortname, allresults=TRUE)
  
   
   
## Base model  
## fit ALL the BASE models without SNPs and save smoking_status-related output to the csv file "Cohort_(GEE/GMMAT)_partial_base_2023.csv"
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq,
                   data_related=data_related,
                   saveOutput=TRUE, cohort=cohortname, BaseModel=TRUE, voi="smoking_status")
  
## fit ALL the BASE models without SNPs and save all the output to the csv file "Cohort_(GEE/GMMAT)_allresults_base_2023.csv"
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq,
                   data_related=data_related,
                   saveOutput=TRUE, cohort=cohortname, BaseModel=TRUE, voi="smoking_status", allresults=TRUE)




## Main effect model for SNP
## Model with only main effect of SNP (no interaction term: SNP*time)

## Save SNP related results:
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq,
                   data_related=data_related,
                   saveOutput=TRUE, cohort=cohortname, SNPmainOnly=TRUE)
   
## Save all results:   
   output <- f_fit(dat=d, covars_common=covars_common, covars_for_eq=covars_for_eq, 
                   data_related=data_related,
                   saveOutput=TRUE, cohort=cohortname, SNPmainOnly=TRUE, allresults=TRUE)















##############################################################  
#            (3). Fit models separately   (Extra)            #
##############################################################
## This section is for analyst to test/debug/explore specific 
## model that does not look correct. 

## Make sure to run Section (1) before running the code below.
## You may either run section (3) or section (2) because 
## section (2) may not work if you run it after running section (3)


## -----------------------------------------------------------
## example code: SNP rs507211

## checking and cleaning the data
 #  d1_test <- check_data(dat=d, covars2=covars_common, snpi="rs507211", summarize=F)   
 #  d2_test <- d_slope(d1_test)                           # Full slope data
 #  d3_test <- d_slope(d1_test,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
 



## Uncomment the code below if you want to test the model manually
## LME
## fit ALL lme models for one SNP:
 #  out_gmmat1 <- fit_all_gmmat(dat_full=d1_test, dat_slope=d2_test, dat_slope_lm=d3_test, eqlist=l_gmmat,          
 #                              covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat)
   
## fit model 3 and model 6 (LME) for one SNP:   
 #  out_gmmat2 <- fit_all_gmmat(dat_full=d1_test, dat_slope=d2_test, dat_slope_lm=d3_test, eqlist=l_gmmat[c(3,6),], 
 #                              covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat) 

## show coefficients for ALL variables
 #  out_gmmat3 <- fit_all_gmmat(dat_full=d1_test, dat_slope=d2_test, dat_slope_lm=d3_test, eqlist=l_gmmat,          
 #                              covars_additional=covars_for_eq, snpi="rs507211", kmat=kmat, all_results=TRUE)


## fit base models for model 1 without SNPs
 #  out_gmmat4 <- fit_all_gmmat(dat_full=d1_test, dat_slope=d2_test, dat_slope_lm=d3_test, eqlist=l_gmmat[1,],      
 #                              covars_additional=covars_for_eq, snpi="smoking_status", kmat=kmat) 
   

 #  out_test_gmmat1
 #  out_test_gmmat2
 #  out_test_gmmat3
 #  out_test_gmmat4



## GEE
## fit ALL GEE models for one SNP
 #  out_gee1 <- fit_all_gee(dat_full=d1_test, dat_slope=d2_test, eqlist=l_gee,           
 #                          covars_additional=covars_for_eq, snpi="rs507211", related=data_related)

## fit model 1 and 3 (GEE) for one SNP
 #  out_gee2 <- fit_all_gee(dat_full=d1_test, dat_slope=d2_test, eqlist=l_gee[c(1,3), ], 
 #                          covars_additional=covars_for_eq, snpi="rs507211", related=data_related)

## show coefficients for ALL variables
 #  out_gee3 <- fit_all_gee(dat_full=d1_test, dat_slope=d2_test, eqlist=l_gee,           
 #                          covars_additional=covars_for_eq, snpi="rs507211", related=data_related, all_results=TRUE)


## fit base models for model 1 without SNPs
 #  out_gee4 <- fit_all_gee(dat_full=d1_test, dat_slope=d2_test, eqlist=l_gee[1, ],      
 #                          covars_additional=covars_for_eq, snpi="smoking_status", related=data_related)


 #  out_test_gee1
 #  out_test_gee2
 #  out_test_gee3
 #  out_test_gee4



