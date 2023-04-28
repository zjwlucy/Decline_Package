# (1) fitting models for multiple SNPs
# (2) write and save the output





##############################################################  
#                 (1) preparation/checking                   #
##############################################################
#  source("f_data.R")
 
#  covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace)
#  covars_for_eq <- covars$covars_for_eq        # use for fitting models
#  covars_common <- covars$covars_common        # use for checking data


# checking and cleaning the data
#  d1 <- check_data(d, covars2=covars_common)   
#  d2 <- d_slope(d1)                           # Full slope data
#  d3 <- d_slope(d1,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
 




##############################################################  
#                       (2) fit models                       #
##############################################################
  library(xlsx)

  f_fit <- function(BaseModel=FALSE, voi=NULL, saveOutput=FALSE, allresults=FALSE){
      
      if(BaseModel){
             l_rs <- voi                                      # for non-SNP variables
      }else{ 
             l_rs <- colnames(d)[grep("^rs", colnames(d))] 
      }
      
      
      all_gmmat <- NULL
      all_gee   <- NULL
      for(tmpSNP in l_rs){
        
        print("#######################")
        print(paste0("SNP:", tmpSNP) )
        print("#######################")
        
        # LME
          tmp_gmmat <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, 
                                     eqlist=l_gmmat, covars_additional=covars_for_eq, snpi=tmpSNP, kmat=kmat, all_results=allresults)
          all_gmmat <- rbind(all_gmmat, tmp_gmmat)
             
        # GEE
          tmp_gee <- fit_all_gee(dat_full=d1, dat_slope=d2, 
                                 eqlist=l_gee, covars_additional=covars_for_eq, snpi=tmpSNP, related=data_related, all_results=allresults)
          all_gee <- rbind(all_gee, tmp_gee)      
      }


     # -------------------    
       all_models  <- list(gmmat=all_gmmat, gee=all_gee) 
       fullresults <- ifelse(allresults, "allresults", "partial")
      
       if(saveOutput){
         if(BaseModel){
               write.xlsx(all_gmmat, file = paste0("summary_",fullresults,"_base_2023.xlsx"), sheetName = "GMMAT", append = TRUE, showNA = F, row.names = F)    
               write.xlsx(all_gee,   file = paste0("summary_",fullresults,"_base_2023.xlsx"), sheetName = "GEE",   append = TRUE, showNA = F, row.names = F)
               save(all_models,      file = paste0("summary_",fullresults,"_base_2023.rdata") )
        
         }else{
               write.xlsx(all_gmmat, file = paste0("summary_",fullresults,"_2023.xlsx"), sheetName = "GMMAT", append = TRUE, showNA = F, row.names = F)    
               write.xlsx(all_gee,   file = paste0("summary_",fullresults,"_2023.xlsx"), sheetName = "GEE",   append = TRUE, showNA = F, row.names = F)
               save(all_models,      file = paste0("summary_",fullresults,"_2023.rdata") )
         }                         
       }
      
     return(all_models) 
     }


