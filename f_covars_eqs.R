## functions to 
## creating covariates list 
## creating equations

## Baseline 
 # age_baseline:           baseline age
 # ht_baseline:            baseline height (in cm)
 # htBaseCenteredSq:       (baseline height - 165)^2
 # smoking_packyears_base: pack-years at baseline
 # sex:                    biological sex
 
## Time-varying
 # smoking_status 
 
## other cohort-specific variables.
 # race
 # PCs 
 # equipchange ......

## Interaction term with time
 # smoking status, sex, smoking_packyears_base, race
      
## grouping variable
 # smoking_status_base: baseline smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin.


## ------------------------------------------------------------------------   
## 
## others: cohort specific covariates

  f_covars <- function(others=NULL,others_inter=NULL,multiRace=FALSE){    
        # common covariates across cohort & will have interaction terms with time
          covars_common   <- c("smoking_status", "sex", "smoking_packyears_base", "ht_baseline", "htBaseCenteredSq") 
          covars_forinter <- c("smoking_status", "sex", "smoking_packyears_base")       # covariates used for interaction with time
                
        # if cohort has race variable     
          if(multiRace){   print("Note: Variable race is added")
                           covars_common   <- c(covars_common,   "race")
                           covars_forinter <- c(covars_forinter, "race")   }
        
        # Interaction terms with time & other interaction terms "others_inter" 
          covars_interWithTime <- c(paste0(covars_forinter, "*", "timefactor_spiro") , others_inter)

    
        # covariates used for fitting models
          covars_for_eq <- c(covars_common, covars_interWithTime, others)
        
        # common covariates used for checking data  
          covars_common <- covars_common[ !(covars_common %in% c("htBaseCenteredSq"))  ]
          covars_common <- c(covars_common, others, "smoking_status_base")
        
        #  
          covars_list   <- list(covars_for_eq, covars_common)
          names(covars_list) <- c("covars_for_eq", "covars_common")
  return(covars_list)     
  }


   
    
      

 
###################################################################### 
###################################################################### 
## load model lists
   library("readxl")
   #l_eq <- read_excel("decline_selectedM_20230109.xlsx")
   l_eq <- read_excel("decline_selectedM_20230424.xlsx")
  
 # replace the name of time^2 by timeCenteredSq
   l_eq$Fixed <- gsub("timefactor_spirosq", "timeCenteredSq", l_eq$Fixed) 
 

## selected models:
 # --------------- GMMAT --------------- 
 # 1	snp*time, age	int,time	snp:time
 # 3  snp*age, time	int,age	snp:age
 # 7	snp*time, age, baseline_ageXtime	int,time	snp:time
 # 9	snp*time, age, timesq	int,time	snp:time
 # 15	snp_s, age	int	snp_s
 # 16 snp_s, age	none	snp_s  (linear model)
   
      # eq1 <- as.formula(paste("fev1", paste(c("snp*time", "age",      "(time|id)"),  collapse=" + "), sep=" ~ "))
      # eq3 <- as.formula(paste("fev1", paste(c("snp*age",  "time",     "(age|id)"),   collapse=" + "), sep=" ~ "))    
      # eq7  <- as.formula(paste("fev1", paste(c("snp*time", "age",   "baseage_time",  "(time|id)"),  collapse=" + "), sep=" ~ "))
      # eq9  <- as.formula(paste("fev1", paste(c("snp*time", "age",      "timesq", "(time|id)"),  collapse=" + "), sep=" ~ "))    
      # eq15 <- as.formula(paste("fev1_s", paste(c("snp_s",    "age",      "(1|id)"), collapse=" + "), sep=" ~ "))
      # eq16 <- as.formula(paste("fev1_s", paste(c("snp_s",    "age"), collapse=" + "), sep=" ~ "))# linear model, outcome = fev1_s
            
   l_gmmat <- l_eq[which(l_eq$modeltype %in% c("glmmkin","lm") & l_eq$model %in% c(1,3,7,9,15,16)), ]
 
 
 # --------------- GEE ---------------
 # 1	 snp*time, age	snp:time
 # 3	 snp*age, time	snp:age
 # 13	 snp_s, age	snp_s
          
      # eq1 <- as.formula(paste("fev1", paste(c("snp*time", "age"),      collapse=" + "), sep=" ~ "))
      # eq3 <- as.formula(paste("fev1", paste(c("snp*age",  "time"),     collapse=" + "), sep=" ~ "))
      # eq13 <- as.formula(paste("fev1_s",    paste(c("snp_s",    "age"),  collapse=" + "), sep=" ~ "))   
      
   l_gee <- l_eq[which(l_eq$modeltype == "GEE" & l_eq$model %in% c(1,3,13)), ]





