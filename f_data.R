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
   check_data <- function(dat, covars2=NULL, snpi=NULL, summarize, cohort_name, rm_change_smk){    
       # basic covariates
         if( !is.null(snpi) && snpi=="smk_status"){ snpi <- "smoking_status" }
         covars1    <- c("IID", "FID", "obsID", "pre_fev1", "timefactor_spiro", "age", "age_baseline", snpi) 
         covars_all <- c(covars1, covars2)   # additional covariates for model (eg: height, race, ......)
         covars_all <- unique(covars_all)
         
       # other variables that are not used in the model
         covars_extra <- base::intersect(c("pre_fev1fvc", "fev1_pp"), colnames(dat)) 
         
       # assign unique id for each observation (used for merging later)  
         dat$obsID <- 1:nrow(dat)
          
       #(1) Missing important covariates
         if( sum(!( covars_all %in% colnames(dat)) ) > 0 ){
                   print("ERROR: missing important covariates: ")
                   print( paste(covars_all[  !(covars_all  %in% colnames(dat))], collapse=",") )
                   #print( paste(covars_extra[!(covars_extra %in% colnames(dat))], collapse=",") )
      
       #(2) Duplicated rows
         }else if( sum(duplicated(dat[,c("IID", "FID", "age")]),
                       duplicated(dat[,c("IID", "FID", "timefactor_spiro")]),
                       duplicated(dat[,covars_all])  
                       ) > 0 ){                                              # check duplicated rows
                   print("ERROR: Data has duplicated rows") 
       
       #(3) Missing values in outcome (FEV1) 
         }else if( sum(is.na(dat$pre_fev1)) >0 ){
                   print("ERROR: Missing values in FEV1") 
         

       #(4) Proceed to next step
         }else{
          # 
            dat <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]
            if( length(grep("\\D", dat$FID)) > 0){ 
                print("WARNING: FID is not numeric, recreating the new FID") 
                n_fid   <- data.frame(FID=unique(dat$FID), new_fid=1:length(unique(dat$FID)) )
                dat     <- merge(dat, n_fid, by="FID", all.x=T) 
                dat$FID <- dat$new_fid
            }
            
            if( length(grep("\\D", dat$IID)) > 0){ 
                print("WARNING: IID is not numeric, recreating the new IID")
                #n_iid   <- as.data.frame(table(dat$IID))
                #dat$IID <- rep(1:length(unique(dat$IID)), n_iid$Freq)    # check identical(dat$IID, rep(unique(dat$IID), n_iid$Freq))
                n_iid   <- data.frame(IID=unique(dat$IID), new_iid=1:length(unique(dat$IID)) )
                dat     <- merge(dat, n_iid, by="IID", all.x=T)
                dat$IID <- dat$new_iid 
            }
  
          #  
            print(paste0("Total number of observations: ", nrow(dat), 
                         "; Number of unique individuals: ", length(unique(dat$IID)) )  )               
            dat$FID                <- as.numeric(dat$FID)
            dat$IID                <- as.numeric(dat$IID)
            dat$smoking_status     <- as.character(dat$smoking_status)
            dat$smoking_status_base<- as.character(dat$smoking_status_base)
            dat$sex                <- as.character(dat$sex)  
          
            dat       <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]
            tmp       <- dat[, c("IID", "obsID", covars_extra)] 
            
          ## ------------- Adding new variable for smoking status 
          ## add a variable indicate consistent smoking status & chaning status  
          ## Note: assuming that missing smoking status will not change the results for consistent status
             const_status <- f_consistent_smk(dat, rm_change_smk_i=FALSE)
             dat          <- merge(dat, const_status, by="IID", all.x=T)                      
             covars_all   <- c(covars_all, "smk_status")
             #table(dat$smoking_status, dat$smk_status)
          ## -------------
          
             
          
          ##----------------------------------------------------------------------  
          ## plots & tables based on dataset WITHOUT genetic information
          ## Using data BEFORE removing missing values  
          ## Assuming that no missing PFTs
             if(summarize){
                print("Generating plots and tables")
                s_table <- f_tables(dat=dat, multiRace=data_multiRace)
                f_plots(dat, cohort_name=cohort_name)
                return(s_table)
                
          ##----------------------------------------------------------------------                  
             }else{
                if(rm_change_smk){ dat <- dat[which(dat$smk_status != "inconsistent"), ] 
                                   print("WARNING: removing individuals with changing smoking status") }
                     
                dat                    <- na.omit(dat[, covars_all]) 
                pft_count              <- as.data.frame(table(dat$IID))      # add number of pfts for each individual
                colnames(pft_count)    <- c("IID", "n_pft")   
                dat                    <- merge(dat, pft_count, by="IID",             all.x=T)
                dat                    <- merge(dat, tmp,       by=c("IID", "obsID"), all.x=T)
              
              # ------------------------ 
                dat$timeCenteredSq     <- (dat$timefactor_spiro - mean(dat$timefactor_spiro, na.rm = T))^2
                dat$htBaseCenteredSq   <- ( dat$ht_baseline - 165 )^2
                #dat$htBaseCenteredSq  <- ( dat$ht_baseline - mean(dat$ht_baseline[which(dat$timefactor_spiro == 0)], na.rm = T) )^2
                #dat$htCenteredSq      <- ( dat$ht_cm       - mean(dat$ht_baseline[which(dat$timefactor_spiro == 0)], na.rm = T) )^2
              # ------------------------
              
              
              # ------------------------ consistent smoking status 
              # # check if smk_status changes after removing missing values
              # !!!!!! smk_status (that are changing) may become consistent after removing missing values.
                const_status           <- f_consistent_smk(dat, rm_change_smk_i=rm_change_smk)
                colnames(const_status) <- c("IID", "smk_status_new")
                dat                    <- merge(dat, const_status, by="IID", all.x=T)                           
                #table(dat$smoking_status, dat$smk_status_new)
                #table(dat$smk_status,     dat$smk_status_new)
                
                if(!identical(dat$smk_status, dat$smk_status_new)){
                   print("WARNING: smoking status changes from inconsistent to consistent or losing 1 level")  
                }
                dat$smk_status <- dat$smk_status_new
              # ------------------------ 
                
                         
                
              # Record the order of each observation for each individuals
                dat_count <- lapply(unique(dat$IID), function(x){
                                    tmp         <- dat[which(dat$IID == x), c("IID", "obsID", "age")]
                                    tmp$obsRank <- rank(tmp$age)
                                    return(tmp)
                                   })
                dat_count <- do.call(rbind, dat_count)
                dat       <- merge(dat, dat_count[,c("IID", "obsID", "obsRank")], by=c("IID", "obsID"), all.x=T)                   
                
                
              # check changes in baseline
                id_base <- which(dat$obsRank == 1)
                check_baseline <- sum(is.logical(  all.equal(as.numeric(dat$age_baseline[id_base]), 
                                                             as.numeric(dat$age[id_base]))  ), 
                                      identical(dat$smoking_status_base[id_base],       dat$smoking_status[id_base])    )                     
                if(check_baseline <2){ print("WARNING: Baseline has changed after removing the missing values") }
         
                
              # ------------------------------------------------ 
                dat <- dat[order(dat$FID, dat$IID, dat$age),]
                print(paste0("Total number of observations after removing NAs: ", nrow(dat), 
                             "; Number of unique individuals: ", length(unique(dat$IID)))  )
                print( paste0("Variables included: ", paste(colnames(dat), collapse=",") )  )                                                  
                
                return(dat)
            }
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
         dat2 <- dat2[order(dat2$FID, dat2$IID, dat2$timefactor_spiro), ]
      
       ## using all observations or not 
         if(firstlast){
             dat2 <- dat2 %>%group_by(IID)%>%slice(c(1,n()))
         }
         
       ## 
         cal_diff <- function(id_i){     
            x   <- dat2[which(dat2$IID == id_i),]
            x   <- x[order(x$FID, x$IID, x$timefactor_spiro), ]
            tmp <- data.frame(x[1:(nrow(x)-1),], 
                              diff(x$pre_fev1), diff(x$timefactor_spiro))
            colnames(tmp) <- c(colnames(x), "fev1_diff", "time_diff")
            return(tmp) 
         }   
            
            
         dat2           <- lapply(unique(dat2$IID), cal_diff)
         dat2           <- do.call("rbind", dat2)    
         dat2$fev1_s    <- dat2$fev1_diff/dat2$time_diff
         rownames(dat2) <- NULL
         colnames(dat2)[grep("^rs", colnames(dat2))] <- paste0(colnames(dat2)[grep("^rs", colnames(dat2))], "_s")
         
         print(paste0("Total number of observations for slope data: ", nrow(dat2), 
                       "; Number of unique individuals: ", length(unique(dat2$IID)))  )
   
         dat2 <- dat2[order(dat2$FID, dat2$IID, dat2$timefactor_spiro), ]
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
## summarize & clean & tools
   f_summary <- function(m_out, forwhich=NULL){
   
           #want  <- c(m_out$variable[grep("^rs", m_out$variable)])  # s$variable[grep(":", s$variable)]
           want  <- m_out$variable[grep(forwhich, m_out$variable)]
           m_out <- m_out[which((m_out$variable %in% want) ), ]
            
   return(m_out) 
   }
   
   
  
## creating smoking variable for consistent status        
   f_consistent_smk <- function(tmp_dat, rm_change_smk_i=FALSE){
         tmp_status <- aggregate(smoking_status~IID, tmp_dat, 
                                 FUN=function(x){ 
                                     x               <- na.omit(x)
                                     changing_status <- length(unique(x))
                                     if(changing_status == 1){      return( unique(x) )
                                     }else if(changing_status > 1){ return("inconsistent") }
                                 })  
       colnames(tmp_status) <- c("IID", "smk_status")                  
       tmp_status$smk_status[which(tmp_status$smk_status == "0")] <- "consistent_none"
       tmp_status$smk_status[which(tmp_status$smk_status == "1")] <- "consistent_former"
       tmp_status$smk_status[which(tmp_status$smk_status == "2")] <- "consistent_current"
       
      #
       if(rm_change_smk_i){ levelist <- c("consistent_none", "consistent_former", "consistent_current") 
       }else{               levelist <- c("consistent_none", "consistent_former", "consistent_current", "inconsistent")} 
       tmp_status$smk_status <- factor(tmp_status$smk_status, levels=levelist)
             
   return(tmp_status)
   }
 
 