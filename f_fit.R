# (1) fitting models for multiple SNPs
# (2) write and save the output





##############################################################  
#                 (1) preparation/checking                   #
##############################################################
#  source("f_data.R")
 
#  covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace)
#  covars_for_eq <- covars$covars_for_eq        # use for fitting models
#  covars_common <- covars$covars_common        # use for checking data





##############################################################  
#                       (2) fit models                       #
##############################################################
  #library(xlsx)
  source("f_plots.R")
  source("f_tables.R")

  # BaseModel   (T/F): TRUE = fit the base model without SNP variable 
  # voi=NULL:          variable of interest. Use only when BaseModel=TRUE
  # saveOutput  (T/F): TRUE = save model results to csv and rdata file
  # allresults  (T/F): TRUE = save all the estimates for all variables (including those that we are not interested in)
  # SNPmainOnly (T/F): TRUE = remove the interaction term for SNP*time and SNP*age, only include main effect of SNP
  f_fit <- function(dat, covars_common, covars_for_eq, data_related,
                    BaseModel=FALSE, voi=NULL, saveOutput=FALSE, allresults=FALSE, cohort="CohortName", SNPmainOnly=FALSE){
      
        if(BaseModel){
              l_rs <- voi                                    # for non-SNP variables
        }else{ 
              l_rs <- colnames(dat)[grep("^rs", colnames(dat))] 
        }
      
      
     ## ------------------------------------------------------------   
     ## plots & tables based on dataset WITHOUT genetic information
     ## Using data BEFORE removing missing values  
     ## Assuming that no missing PFTs
        s_table <- check_data(dat=dat, covars2=covars_common, snpi=NULL, summarize=TRUE, cohort_name=cohort)  

          
     ## ------------------------------------------------------------  
     ## Fit model
        all_gmmat <- NULL
        all_gee   <- NULL
        for(tmpSNP in l_rs){
        
          print("#############################")
          print(paste0("SNP:", tmpSNP) )
          print("#############################")
        
          # checking and cleaning the data
            d1 <- check_data(dat=dat, covars2=covars_common, snpi=tmpSNP, summarize=FALSE)   
            d2 <- d_slope(d1)                           # Full slope data
            d3 <- d_slope(d1,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
 
 
          # LME
            tmp_gmmat <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, snpi=tmpSNP, SNP_mainOnly=SNPmainOnly,
                                       eqlist=l_gmmat, covars_additional=covars_for_eq, kmat=kmat, all_results=allresults)
            all_gmmat <- rbind(all_gmmat, tmp_gmmat)
             
          # GEE
            tmp_gee <- fit_all_gee(dat_full=d1, dat_slope=d2, snpi=tmpSNP, SNP_mainOnly=SNPmainOnly,
                                   eqlist=l_gee, covars_additional=covars_for_eq, related=data_related, all_results=allresults)
            all_gee <- rbind(all_gee, tmp_gee)      
        }

              
      
     ## -------------------   
     ## save output 
        all_models  <- list(s_table=s_table, gmmat=all_gmmat, gee=all_gee) 
        fullresults <- ifelse(allresults, "allresults", "partial")
      
        if(saveOutput){
          if(BaseModel){
                write.csv(all_gmmat, file=paste0(cohort, "_GMMAT_",  fullresults,"_base_2023.csv"), row.names = F)   
                write.csv(all_gee,   file=paste0(cohort, "_GEE_",    fullresults,"_base_2023.csv"), row.names = F)
                save(all_models,     file=paste0(cohort, "_summary_",fullresults,"_base_2023.rdata") )
        
          }else{
                SNPeffect <- ifelse(SNPmainOnly, "SNPmain_", "")
                write.csv(all_gmmat, file=paste0(cohort, "_GMMAT_",  fullresults,"_", SNPeffect, "2023.csv"), row.names = F)    
                write.csv(all_gee,   file=paste0(cohort, "_GEE_",    fullresults,"_", SNPeffect, "2023.csv"), row.names = F)
                save(all_models,     file=paste0(cohort, "_summary_",fullresults,"_", SNPeffect, "2023.rdata") )
          }
          write.csv(s_table,   file=paste0(cohort, "_table_2023.csv"), row.names = F)                           
        }
      
      

     return(all_models) 
     }







