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
          
           # related individuals (cluster on FID)
             if(relatedi){
                 m_id   <- "FID"
                 dat    <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]
                 timei  <- system.time({  m_gee <- geeglm(eqi, id=FID, corstr="unstructured", data=dat)  })          
             
           # unrelated individuals (cluster on IID)
             }else{
                 m_id   <- "IID"
                 dat    <- dat[order(dat$IID, dat$timefactor_spiro),]
                 timei  <- system.time({  m_gee <- geeglm(eqi, id=IID, corstr="unstructured", data=dat)  })    
             }
             
           # 
             s <- as.data.frame(summary(m_gee)$coefficients)
             s <- cbind(modeltypei, modeli, dim(dat)[1], length(unique(dat$IID)), rownames(s), s, 
                        QIC(m_gee)[1], as.numeric(timei["sys.self"]), as.numeric(timei["elapsed"]), m_id, relatedi)
             colnames(s) <- c("modeltype", "model", "n_obs", "n_uniq", "variable", "Estimate", "SE", "Wald", "pvalue",
                              "QIC", "sys_time", "elapsed_time", "m_id","related")

   return(s)
   }
          
          
          
          
     
## ----------------------------------------------------------
## dat:    full data
## dat2:   full slope data

## eqlist:      list of models to fit
## related:     does the data contain related or unrelated individuals
## all_results: whether to present coefficient for all 

   fit_all_gee <- function(dat_full, dat_slope, eqlist, covars_additional=NULL, snpi, related, all_results=FALSE, SNP_mainOnly=FALSE){
       
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
             # SNP main effect model or not
               if(SNP_mainOnly){ covarsi <- gsub("snp[*]timefactor_spiro", "snp, timefactor_spiro", covarsi)
                                 covarsi <- gsub("snp[*]age",              "snp, age",              covarsi) }
                         
             # construct formula for specified SNP or variable
               if(eqlist$variable[i] == "snp_s"  &  !grepl("^rs", snpi) ){
                    covarsi <- gsub("snp_s", snpi, covarsi)     # for non-SNP variables
               }else{
                    covarsi <- gsub("snp", snpi, covarsi)       # for SNPs 
               }
             # Remove the interaction with time for slope data: 
               if(eqlist$variable[i] == "snp_s"){
                  covars_additional <- gsub("[*]timefactor_spiro", "", covars_additional)
               }
             #     
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

      
 
          
      

       