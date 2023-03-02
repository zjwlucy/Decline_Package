##############################################################  
#                  GMMAT (glmmkin/lmekin)                    #
##############################################################  
# ------------------------------------------------------------
## selected models:
## GMMAT: 
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
        
    
  #l_gmmat <- l_eq[which(l_eq$modeltype == "glmmkin" & l_eq$model %in% c(1,3,7,9,15,16)), ]
 



##############################################################  
# ------------------------------------------------------------
# dat: dataset
# modeli: ith model
# IID: individual id

  library(GMMAT)
  f_glmmkin <- function(dat, modeli, eqi, rand_s=NULL, m_id="IID", m_group="smoking_status", kmati=NULL, modeltypei=NULL){
        use_kmat   <- FALSE
        if(!is.null(kmati)){  
          ids      <- rownames(kmati)[which(rownames(kmati) %in% dat$IID)]  
          indexi   <- !is.na(match(dimnames(kmati)[[1]], ids))    
          kmati    <- kmati[indexi, indexi]
          use_kmat <- TRUE
        }
      #      
        timei  <- system.time({ m_gmmat <- glmmkin(eqi, random.slope = rand_s, 
                                                   id = m_id, groups = m_group, kins = kmati, 
                                                   data = dat, family = gaussian(link = "identity"))   }) 
      #                                                                              
        fixed   <- as.data.frame( cbind( m_gmmat$coefficients, sqrt(diag(m_gmmat$cov)) ) )   
        fixed$t <- fixed$V1/fixed$V2
        fixed$p <- 2*pnorm(abs(fixed$V1/fixed$V2), lower.tail = F)       
        colnames(fixed) <- c("Estimate", "SE", "t", "pvalue") 
                
  
      #
        if(is.null(m_group)){ m_group <- "NONE" }
        s           <- cbind(modeltypei, modeli, dim(dat)[1], length(unique(dat$IID)), rownames(fixed), fixed, timei[2], timei[3], use_kmat, m_id, m_group)   
        colnames(s) <- c("modeltype", "model", "n_obs", "n_uniq", "variable", "Estimate", "SE", "t", "pvalue", 
                         "sys_time", "elapsed_time", "UseKinship", "m_id", "m_group") 
           
   return(s)             
   }
 



 
## ------------------------------------------------------------
## dat:  full data
## dat2: full slope data
## dat3: slope data with 1 observation

## eqlist:      list of models to fit
## all_results: whether to present coefficient for all 

   fit_all_gmmat <- function(dat_full, dat_slope, dat_slope_lm, eqlist, 
                             covars_additional=NULL, snpi, kmat=NULL, all_results=FALSE){
   
     # ----------------------------     
       gmmat_out <- NULL 
       for(i in 1:nrow(eqlist) ){ 
         tryCatch({      
             modeltypei <- eqlist$modeltype[i]
             modeli     <- eqlist$model[i]
             covarsi    <- eqlist$Fixed[i]
             randi      <- eqlist$Rand[i]
             outcomei   <- ifelse(eqlist$variable[i] == "snp_s", "fev1_s", "pre_fev1")
             if(randi == "int" | is.na(randi) ){ randi <- NULL}
             print("=========================================================================")
             print( paste0("Fitting ", modeltypei, " model ", modeli, " with random slope (", randi, ") and outcome ", outcomei) )
             
           # construct formula for specified SNP or variable
             if(eqlist$variable[i] == "snp_s"  &  !grepl("^rs", snpi) ){
                   covarsi <- gsub("snp_s", snpi, covarsi)     # for non-SNP variables
             }else{
                   covarsi <- gsub("snp", snpi, covarsi)       # for SNPs 
             }
             covarsi <- gsub("timefactor_spirosq", "timeCenteredSq", covarsi) 
             covarsi <- c(unlist(strsplit(covarsi, split = ", ")), covars_additional) 
             covarsi <- unique(covarsi)
             covarsi <- paste(covarsi, collapse=" + ")       
             eqi     <- as.formula(paste(outcomei, covarsi, sep=" ~ ")) 
             print( paste0("Fixed terms: ", covarsi ) )
            
           # -------------          
           # (A) linear model using slope data (with only 1 observation per individual, no random slope)
           #    (1) 1 observation, independent individuals  (ID=IID, no kmat, no slope)   we cannot use m_group="smoking_status"
           #    (2) 1 observation, related individuals  (ID=IID, kmat, no slope)    ????do we need m_group="smoking_status"
                if(eqlist$variable[i] == "snp_s" & modeltypei == "lm"){  
                  
                   if(is.null(kmat)){ m_groupi <- NULL }else{ m_groupi <- "smoking_status" }            
                   m_tmp <- f_glmmkin(dat_slope_lm, modeli, eqi, rand_s=NULL, m_group=m_groupi, kmati=kmat, modeltypei=modeltypei)
               
           # (B) lme model using slope data (multiple slopes as outcome, no random slope)  
                 }else if(eqlist$variable[i] == "snp_s" & modeltypei == "glmmkin"){              
                  m_tmp <- f_glmmkin(dat_slope, modeli, eqi, rand_s=NULL, kmati=kmat, modeltypei=modeltypei)
             
           # (C) lme models    
                 }else{     
                  m_tmp <- f_glmmkin(dat_full, modeli, eqi, rand_s=randi, kmati=kmat, modeltypei=modeltypei)
                 }
           # -------------
             
             if(is.null(randi)){ randi <- "NONE"}
             m_tmp         <- cbind(m_tmp, SNP=snpi, randslope=randi, outcome=outcomei)
             m_tmp$covarsi <- covarsi
             m_tmp$coefs   <- paste(m_tmp$variable, collapse=",")  
                                                                    
             
             gmmat_out <- rbind(gmmat_out,  m_tmp) 
        
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
       }
     # ----------------------------  
   
     # whether to save coefficients for all variables
        if(!all_results){ 
                gmmat_out <- f_summary(gmmat_out, forwhich = snpi)
        } 
       gmmat_out <-  gmmat_out[order(gmmat_out$SNP, gmmat_out$modeltype, gmmat_out$model), ]
       rownames(gmmat_out) <- NULL
 
   return(gmmat_out)
   }






# linear model with unrelated individuals
#  m_gmmat <- glmmkin(eqi, random.slope = NULL, id = "IID", groups = NULL, kins = NULL,
#                     data = dat_slope_lm, family = gaussian(link = "identity"))

# linear model with unrelated individuals + hetero  (does not work)                   
#  m_gmmat2 <- glmmkin(eqi, random.slope = NULL, id = "IID", groups = "smoking_status", kins = NULL,
#                     data = dat_slope_lm, family = gaussian(link = "identity"))
                     

# linear model with related individuals + hetero    
# use "IID" or "FID"????                    
#  m_gmmat2 <- glmmkin(eqi, random.slope = NULL, id = "IID", groups = "smoking_status", kins = kmati,
#                     data = dat_slope_lm, family = gaussian(link = "identity"))                   
                     


