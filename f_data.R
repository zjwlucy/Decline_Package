## functions to 
## (1) roughly checking dataset
## (2) generate specific data for models

        
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## covariates
 
## Must have: IID, FID, pre_fev1, SNPs, timefactor_spiro, age, smoking_status (never, former, current)

## Baseline: age_baseline, ht_baseline (cm), htBaseCenteredSq (baseline height - 165)^2,
 # smoking_packyears_base, sex, race, PCs 
 # other cohort-specific variables.

## Interaction term with time
 #   ht_baseline , htBaseCenteredSq, race, sex, smoking_packyears_base, smoking status


## 
   check_data <- function(dat, covars2=NULL){    
       # basic 
         covars1 <- c("IID", "FID", "pre_fev1", colnames(dat)[grep("^rs", colnames(dat))], 
                      "timefactor_spiro", "age", "age_baseline") 
       # additional covariates
         covars_all <- c(covars1, covars2)   
       
       #  
         if( sum(!(covars_all %in% colnames(dat))) > 0 ){
            print("missing covariates: ")
            print( paste(covars_all[!(covars_all %in% colnames(dat))], collapse=",") )
             
         }else{
            print(paste0("Total number of observations: ", nrow(dat), 
                         "; Number of unique individuals: ", length(unique(dat$IID)))  )
            
            dat                    <- na.omit(dat[, covars_all]) 
            pft_count              <- data.frame(table(dat$IID))      # add number of pfts for each individual
            colnames(pft_count)    <- c("IID", "n_pft")   
            dat                    <- merge(dat, pft_count, by = "IID", all.x=T)
          
          # ------------------------ 
            dat$smoking_status     <- as.character(dat$smoking_status)
            dat$sex                <- as.character(dat$sex)  
            dat$timeCenteredSq     <- (dat$timefactor_spiro - mean(dat$timefactor_spiro, na.rm = T))^2
            #dat$htCenteredSq     <- ( dat$ht_cm       - mean(dat$ht_baseline[which(dat$timefactor_spiro == 0)], na.rm = T) )^2
            #dat$htBaseCenteredSq <- ( dat$ht_baseline - mean(dat$ht_baseline[which(dat$timefactor_spiro == 0)], na.rm = T) )^2
            dat$htBaseCenteredSq   <- ( dat$ht_baseline - 165 )^2
          # ------------------------
            
            dat <- dat[order(dat$FID, dat$IID, dat$age),]
            print(paste0("Total number of observations after removing NAs: ", nrow(dat), 
                         "; Number of unique individuals: ", length(unique(dat$IID)))  )
            print( paste0("Variables included: ", paste(colnames(dat), collapse=",") )  )                                                  
            
            return(dat)
         }
        }




## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## generate slope data  
## (1) using all observations
## (2) using first-last observations 
   library(dplyr)
   d_slope <- function(dat, firstlast=FALSE){
         print( paste0("Number of individuals with 1 observation is ", sum(dat$n_pft < 2)) )
         dat2 <- dat[which(dat$n_pft > 1), ]
         dat2 <- dat2[order(dat2$IID, dat2$n_pft, dat2$age), ]
      
       # using all observations or not 
         if(firstlast){
             dat2 <- dat2 %>%group_by(IID)%>%slice(c(1,n()))
         }
       # 
         cal_diff <- function(id_i){     
            x   <- dat2[which(dat2$IID == id_i),]
            tmp <- data.frame(x[1:(nrow(x)-1),], 
                              diff(x$pre_fev1), diff(x$timefactor_spiro))
            colnames(tmp) <- c(colnames(x), "fev1_diff", "time_diff")
            colnames(tmp)[]
            return(tmp) 
         }   
            
            
         dat2           <- lapply(unique(dat2$IID), cal_diff)
         dat2           <- do.call("rbind", dat2)    
         dat2$fev1_s    <- dat2$fev1_diff/dat2$time_diff
         rownames(dat2) <- NULL
         colnames(dat2)[grep("^rs", colnames(dat2))] <- paste0(colnames(dat2)[grep("^rs", colnames(dat2))], "_s")
         
         print(paste0("Total number of observations for slope data: ", nrow(dat2), 
                       "; Number of unique individuals: ", length(unique(dat2$IID)))  )
   
         dat2 <- dat2[order(dat2$FID, dat2$IID, dat2$age), ]
         return(dat2)
        } 

 

  
## first and last observations only   
#    d_2bs <- function(dat){
#           d_tmp <- dat[which(dat$n_pft > 1), ]
#           d_tmp <- d_tmp %>%group_by(IID)%>%slice(c(1,n()))
#    return(d_tmp)
#    }
  
 
 
 
 
## ------------------------------------------------------------------- 
## 
   source("f_gee.R")
   source("f_gmmat.R")
   source("f_covars_eqs.R")









###################################################################### 
###################################################################### 
## summarize & clean
   f_summary <- function(m_out){
   
           want  <- c(m_out$variable[grep("^rs", m_out$variable)])  # s$variable[grep(":", s$variable)]
           m_out <- m_out[which((m_out$variable %in% want) ), ]
            
   return(m_out) 
   }
        
 
 
 