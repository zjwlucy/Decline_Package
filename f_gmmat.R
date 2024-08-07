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
  f_glmmkin <- function(dat, modeli, eqi, rand_s=NULL, m_id="IID", m_group, kmati=NULL, modeltypei=NULL){
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
                
      # calculate MAF
        haveSNP <- grepl("^rs", colnames(m_gmmat$X))
        if( sum(haveSNP)>0 ){
        
           haveSNP <- colnames(m_gmmat$X)[haveSNP]
           haveSNP <- haveSNP[!grepl(":", haveSNP)]   # remove interaction terms
           mafi    <- lapply(unique(m_gmmat$id_include), function(x){ unique(dat[which(dat$IID == x), haveSNP]) } )
           
           if( sum(unlist(  lapply(mafi, function(x){length(x)>1})  )) >0 ){
               stop("ERROR: SNP variable is time varying")
           }else{
               mafi <- unlist(mafi)
               mafi <- sum(mafi)/(2*length(mafi))
               mafi <- min(mafi, 1-mafi)
           } 
        }else{ 
           mafi <- NA
        }         
                
      # get estimates for variance components
      # For lm:
      # kmat group
      #  T    T
      #  F    F
      
      # For LME
      # kmat group
      #  T    T
      #  F    T
        if(modeltypei == "lm" && is.null(kmati)){  names(m_gmmat$theta) <- "1" }  # no random & group specific variance estimate
        vars_info <- data.frame(matrix(vector(), 12, 0,
                                dimnames=list(c("1", "2", "3", "4", 
                                              "kins1.var.intercept", "kins2.var.intercept", 
                                              "kins1.var.slope", "kins2.var.slope",
                                              "kins1.cov.intercept.slope", "kins2.cov.intercept.slope",
                                              "kins1", "kins2"), c()) ),
                                stringsAsFactors=F)
        vars_info$var_name <- rownames(vars_info)
        vars_i             <- data.frame(m_gmmat$theta)
        vars_i$var_name    <- rownames(vars_i)
        vars_info          <- merge(vars_info, vars_i, by = "var_name", all.x=T) 
        vars_info$var_name <- paste0("V_", vars_info$var_name)
        
      #        
        if(is.null(m_group)){  m_group <- "NONE"  }        
        s           <- cbind(modeltypei, modeli, length(m_gmmat$Y), length(unique(m_gmmat$id_include)), rownames(fixed), fixed, 
                             mafi, as.numeric(timei["sys.self"]), as.numeric(timei["elapsed"]), use_kmat, m_id, m_group, 
                             data.frame(t(matrix(vars_info$m_gmmat.theta, nrow=nrow(vars_info), ncol=nrow(fixed))))  )   
        colnames(s) <- c("modeltype", "model", "n_obs", "n_uniq", "variable", "Estimate", "SE", "t", "pvalue", 
                         "MAF", "sys_time", "elapsed_time", "UseKinship", "m_id", "m_group", vars_info$var_name) 
           
   return(s)             
   }
 



 
## ------------------------------------------------------------
## dat:  full data
## dat2: full slope data
## dat3: slope data with 1 observation

## eqlist:      list of models to fit
## all_results: whether to present coefficient for all 

   fit_all_gmmat <- function(dat_full, dat_slope, dat_slope_lm, eqlist, snpi, covars_additional=NULL,
                             m_groupi, kmat=NULL, SNP_mainOnly=FALSE, all_results=FALSE){
   
     # -----------------------------------     
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
             
           ## 
           # (1) SNP main effect model or not
             if(SNP_mainOnly){ covarsi <- gsub("snp[*]timefactor_spiro", "snp, timefactor_spiro", covarsi)
                               covarsi <- gsub("snp[*]age",              "snp, age",              covarsi) }
                         
           # (2) construct formula for specified SNP or variable
             if(!grepl("^rs", snpi) && eqlist$variable[i] == "snp_s"){
                                  covarsi <- gsub("snp_s",  snpi, covarsi)    # for non-SNP variables (base model) in slope models                  
             }else if(!grepl("^rs", snpi) && modeli==3 ){  
                                  covarsi <- gsub("snp[*]", "",   covarsi)    # for non-SNP variable (base model) from model 3
             }else{               covarsi <- gsub("snp",    snpi, covarsi)    # for SNPs 
             }
               
           # (3) Remove interaction terms with time for slope data: 
             if(eqlist$variable[i] == "snp_s"){
                 covars_additional <- gsub("[*]timefactor_spiro", "", covars_additional)
             }
             
           # (4)  
             covarsi <- c(unlist(strsplit(covarsi, split = ", ")), covars_additional) 
             covarsi <- unique(covarsi)
             covarsi <- paste(covarsi, collapse=" + ")       
             eqi     <- as.formula(paste(outcomei, covarsi, sep=" ~ ")) 
             print( paste0("Fixed terms: ", covarsi ) )
             
            
           ## ------------------------          
           # (A) linear model using slope data (with only 1 observation per individual, no random slope)
           #    (1) 1 observation, independent individuals  (ID=IID, no kmat, no slope)   we cannot use m_group="smoking_status_base" when no kmat
           #    (2) 1 observation, related individuals      (ID=IID, kmat,    no slope)   ????do we need m_group="smoking_status_base"? Default is yes.
                if(eqlist$variable[i] == "snp_s" && modeltypei == "lm"){  
                  
                   if(is.null(kmat)){ m_groupi <- NULL }#else{ m_groupi <- "smoking_status_base" }            
                   m_tmp <- f_glmmkin(dat_slope_lm, modeli, eqi, rand_s=NULL, m_group=m_groupi, kmati=kmat, modeltypei=modeltypei)
            
            
           # For LME, m_group will always be m_groupi (baseline, time-varying or consistent smoking status)    
           # (B) lme model using slope data (multiple slopes as outcome, no random slope)  
                 }else if(eqlist$variable[i] == "snp_s" && modeltypei == "glmmkin"){       
                # check if the slope model is fitted on unrelated individuals with 1 observation
                  if(is.null(kmat)){
                       #n_size        <- as.data.frame( table(dat_slope$IID) ) 
                       #change_groupi <- (sum(n_size$Freq==1) == nrow(n_size)) 
                       change_groupi <- ( length(unique(dat_slope$IID)) == nrow(dat_slope) ) 
                       if(change_groupi){
                          m_groupi <- NULL 
                          print("WARNING: dataset is not longitudinal, all individuals have only 1 observation")
                       }
                  }      
                  m_tmp <- f_glmmkin(dat_slope, modeli, eqi, rand_s=NULL, m_group=m_groupi, kmati=kmat, modeltypei=modeltypei)
             
             
           # (C) lme models    
                 }else{     
                  m_tmp <- f_glmmkin(dat_full, modeli, eqi, rand_s=randi, m_group=m_groupi, kmati=kmat, modeltypei=modeltypei)
                 }
           ## ------------------------       
             
             if(is.null(randi)){ randi <- "NONE"}
             m_tmp         <- cbind(m_tmp, SNP=snpi, randslope=randi, outcome=outcomei)
             m_tmp$covarsi <- covarsi
             m_tmp$coefs   <- paste(m_tmp$variable, collapse=",")  
                                                                    
             
             gmmat_out <- rbind(gmmat_out,  m_tmp) 
        
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
       }
     # -----------------------------------  
   
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
                     


