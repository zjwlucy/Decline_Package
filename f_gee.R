############################################################## 
#                            GEE                             #
##############################################################  
# ------------------------------------------------------------
## selected models:
## GEE:
 #1	  snp*time, age	snp:time
 #3	  snp*age, time	snp:age
 #13	snp_s, age	snp_s
          
      # eq1 <- as.formula(paste("fev1", paste(c("snp*time", "age"),      collapse=" + "), sep=" ~ "))
      # eq3 <- as.formula(paste("fev1", paste(c("snp*age",  "time"),     collapse=" + "), sep=" ~ "))
      # eq13 <- as.formula(paste("fev1_s",    paste(c("snp_s",    "age"),  collapse=" + "), sep=" ~ "))   
      
  #l_gee <- l_eq[which(l_eq$modeltype == "GEE" & l_eq$model %in% c(1,3,13)), ]

 

##############################################################  
## ------------------------------------------------------------
## dat: dataset
## modeli: ith model
## relatedi: does the data contain related or unrelated individuals

   library(geepack)
   f_gee <- function(dat, modeli, eqi, relatedi, modeltypei=NULL){    
                     
           # if related, calculate the #obs by FID
             if(relatedi){   n_size <- as.data.frame( table(dat$FID) )
             }else{          n_size <- as.data.frame( table(dat$IID) )     }
             
             
           # check if the dataset is still a longitudinal data 
             alleq1     <- sum(n_size$Freq == 1) == nrow(n_size)
             if(alleq1) { print("WARNING: dataset is not longitudinal, all individuals have only 1 observation") }
             
           # If all individuals have (#obs <= 3) or 90% of them have (#obs <= 3),
           # change unstructured to ar1
             change_cor <- ( sum(n_size$Freq <= 3) == nrow(n_size) | sum(n_size$Freq <= 3) >= 0.90*nrow(n_size)  )
             
             if(change_cor){  cor_str <- "ar1";  print("WARNING: observations per individual is <= 3, changing corstr to AR1 for GEE")
             }else{           cor_str <- "unstructured"      }
            
            
           # check if smoking status have levels with 0 observation (GEE cannot run when factors have used levels)
             if(sum(table(dat$smk_status) > 0) ){ dat$smk_status <- factor(dat$smk_status) }
            
          
           # -----------------------------------------------
           # related individuals (cluster on FID)
             if(relatedi){
                 m_id   <- "FID"
                 dat    <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]
                 timei  <- system.time({  m_gee <- geeglm(eqi, id=FID, corstr=cor_str, data=dat)  })          
             
           # unrelated individuals (cluster on IID)
             }else{
                 m_id   <- "IID"
                 dat    <- dat[order(dat$IID, dat$timefactor_spiro),]
                 timei  <- system.time({  m_gee <- geeglm(eqi, id=IID, corstr=cor_str, data=dat)  })    
             }
             
           # calculate MAF 
              haveSNP <- grepl("^rs", colnames(m_gee$geese$X))
              if( sum(haveSNP)>0 ){
        
                 haveSNP <- colnames(m_gee$geese$X)[haveSNP]
                 haveSNP <- haveSNP[!grepl(":", haveSNP)]   # remove interaction terms
                 mafi    <- lapply(unique(m_gee$geese$id), function(x){ unique(dat[which(dat$IID == x), haveSNP]) } )
                 
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
              
           # 
             s <- as.data.frame(summary(m_gee)$coefficients)
             s <- cbind(modeltypei, modeli, dim(dat)[1], length(unique(dat$IID)), rownames(s), s, 
                        mafi, QIC(m_gee)[1], as.numeric(timei["sys.self"]), as.numeric(timei["elapsed"]), m_id, relatedi)
             colnames(s) <- c("modeltype", "model", "n_obs", "n_uniq", "variable", "Estimate", "SE", "Wald", "pvalue",
                              "MAF", "QIC", "sys_time", "elapsed_time", "m_id","related")

   return(s)
   }
          
          
          
          
     
## ----------------------------------------------------------
## dat:    full data
## dat2:   full slope data

## eqlist:      list of models to fit
## related:     does the data contain related or unrelated individuals
## all_results: whether to present coefficient for all 

   fit_all_gee <- function(dat_full, dat_slope, snpi, eqlist, covars_additional=NULL, related, SNP_mainOnly=FALSE, all_results=FALSE){
       
       # ----------------------------  
         gee_out <- NULL 
         for(i in 1:nrow(eqlist) ){ 
           tryCatch({      
               modeltypei <- eqlist$modeltype[i]
               modeli     <- eqlist$model[i]
               covarsi    <- eqlist$Fixed[i]
               outcomei   <- ifelse(eqlist$variable[i] == "snp_s", "fev1_s", "pre_fev1")
               print("=========================================================================")
               print( paste0("Fitting ", modeltypei, " model ", modeli) )
           
             ##  
             # (1) SNP main effect model or not
               if(SNP_mainOnly){ covarsi <- gsub("snp[*]timefactor_spiro", "snp, timefactor_spiro", covarsi)
                                 covarsi <- gsub("snp[*]age",              "snp, age",              covarsi) }
                         
             # (2) construct formula for specified SNP or variable
               if(!grepl("^rs", snpi) && eqlist$variable[i] == "snp_s"){
                                   covarsi <- gsub("snp_s",  snpi, covarsi)    # for non-SNP variables (base model) in slope models                  
               }else if(!grepl("^rs", snpi) && modeli==3 ){  
                                   covarsi <- gsub("snp[*]", "",   covarsi)    # for non-SNP variable (base model) from model 3
               }else{              covarsi <- gsub("snp",    snpi, covarsi)    # for SNPs 
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
              
             ##--------------------------------              
             # GEE using slope data
               if(eqlist$variable[i] == "snp_s"){              
                  m_tmp <- f_gee(dat_slope, modeli, eqi, relatedi=related, modeltypei=modeltypei)           
             # GEE using full data
               }else{
                  m_tmp <- f_gee(dat_full, modeli, eqi, relatedi=related, modeltypei=modeltypei)
               }  
             #    
               m_tmp         <- cbind(m_tmp, SNP=snpi, outcome=outcomei)
               m_tmp$covarsi <- covarsi
               m_tmp$coefs   <- paste(m_tmp$variable, collapse=",")  
                                       
               gee_out <- rbind(gee_out,  m_tmp) 
               
           }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
         }    
      ## ----------------------------          
      
      # whether to save coefficients for all variables
        if(!all_results){ 
                gee_out <- f_summary(gee_out,forwhich = snpi)
        } 
        gee_out <-  gee_out[order(gee_out$SNP, gee_out$modeltype, gee_out$model), ]
        rownames(gee_out) <- NULL
    
   return(gee_out)
   }

      
 
          
      

       