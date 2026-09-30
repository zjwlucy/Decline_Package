# (1) fitting models for multiple SNPs
# (2) write and save the output



##############################################################  
#                       (2) fit models                       #
##############################################################
  #library(xlsx)
  source("f_data.R")
  source("f_plots.R")
  source("f_tables.R")

  # saveOutput  (T/F):   TRUE = save model results to csv and rdata file
  # allresults  (T/F):   TRUE = save all the estimates for all variables (including those that we are not interested in)
  
  # BaseModel   (T/F):   TRUE = fit the base model without SNP variable 
  # voi=NULL:            variable of interest. Use only when BaseModel=TRUE
  # SNPmainOnly (T/F):   TRUE = remove the interaction term for SNP*time and SNP*age, only include main effect of SNP
  # fixed_smk   (T/F):   TRUE = use newly created smoking variable (smk_status, indicating consistent smoking status) as the smoking variable
  # rm_change_smk (T/F): TRUE = removing individuals with changing smoking status
  
  f_fit <- function(dat, covars_cohort, data_related, cohort="CohortName", saveOutput=TRUE, allresults=FALSE, 
                    BaseModel=FALSE, voi=NULL, SNPmainOnly=FALSE, fixed_smk, rm_change_smk=FALSE){
       
      # create output folder 
        outfolder <- "decline_package_output"
        if(!dir.exists(outfolder)){
            print("Creating output folder") 
            dir.create(outfolder, recursive = TRUE)
        }
      
      
      #----------------------------------------------
      # Create covariate list
      #----------------------------------------------
        covars        <- f_covars(others=covars_cohort,multiRace=data_multiRace,fixed_smk=fixed_smk)
        covars_for_eq <- covars$covars_for_eq        # use for fitting models
        covars_common <- covars$covars_common        # use for checking data
        rs_want       <- covars$rs_want              # list of selected SNPs/variables
        print("Selected SNPs: "); print(rs_want) 
      
      
      # ------------------
      # Main variables: 
        l_rs <- colnames(dat)[grep("^rs", colnames(dat))] 
        
      # 1) Base model without SNPs (for non-SNP variables)
        if(BaseModel){
             l_rs <- voi                              
  
      # 2) Model with SNPs
        }else{
        
         # Model with only SNP main effect
           if(SNPmainOnly){
              l_gmmat <- l_gmmat[which(l_gmmat$variable != "snp_s"), ]
              l_gee   <- l_gee[  which(l_gee$variable   != "snp_s"), ]
           }  
              
         # Model with SNPxTime  
           if(length(l_rs) == 0){                          # Missing SNP columns
              stop("ERROR: missing all pre-selected SNPs") 
         
           }else if( sum(rs_want %in% l_rs) > 0 ){         # Have SNPs that are pre-selected
              
              if( sum(rs_want %in% l_rs) < length(rs_want) ){
                      print("WARNING: missing some SNP columns; Fitting models on SNPs that are presented in the data")
                      l_rs <- base::intersect(l_rs, rs_want)
              }else{  l_rs <- rs_want  } 
             
           }else{                                          # other scenario
              stop("ERROR: check SNP columns")
           }
        }

      # smoking variables for grouping (heterogeneous variances)
        #m_groupi <- ifelse(fixed_smk==TRUE, "smk_status", "smoking_status_base")
         m_groupi <- "smk_status"
     
      
     ## ------------------------------------------------------------   
     ## plots & tables based on dataset WITHOUT genetic information
     ## Using data BEFORE removing missing values  
     ## Assuming that no missing PFTs
        s_table <- check_data(dat=dat, covars2=covars_common, snpi=NULL, summarize=TRUE, cohort_name=cohort, rm_change_smk=FALSE)  

          
     ## ------------------------------------------------------------  
     ## Fit model
        all_gmmat <- NULL
        all_gee   <- NULL
        for(tmpSNP in l_rs){
        
          print("#############################")
          print(paste0("SNP:", tmpSNP) )
          print("#############################")
        
          # checking and cleaning the data
            d1 <- check_data(dat=dat, covars2=covars_common, snpi=tmpSNP, summarize=FALSE, rm_change_smk=rm_change_smk)   
            d2 <- d_slope(d1)                           # Full slope data
            d3 <- d_slope(d1,firstlast=T)               # single slope data where slope = (last-first)/(time interval)
          # table(d1$smoking_status, d1$smk_status)
 
 
          # LME
            tmp_gmmat <- fit_all_gmmat(dat_full=d1, dat_slope=d2, dat_slope_lm=d3, snpi=tmpSNP, eqlist=l_gmmat, 
                                       covars_additional=covars_for_eq, m_groupi=m_groupi, kmat=kmat,
                                       SNP_mainOnly=SNPmainOnly, all_results=allresults)
            all_gmmat <- rbind(all_gmmat, tmp_gmmat)
             
             
          # GEE
            tmp_gee   <- fit_all_gee(dat_full=d1, dat_slope=d2, snpi=tmpSNP, eqlist=l_gee,
                                     covars_additional=covars_for_eq, related=data_related, 
                                     SNP_mainOnly=SNPmainOnly, all_results=allresults)
            all_gee   <- rbind(all_gee, tmp_gee)      
        }
             
      
     ## -------------------   
     ## save output 
        all_models  <- list(s_table=s_table, gmmat=all_gmmat, gee=all_gee) 
        
      #            
        if(saveOutput){
        
          if(fixed_smk){       fixed_smk <- "newSMK_"             # consistent + inconsistent
            if(rm_change_smk){ fixed_smk <- "exInconsistSMK_"  }  # only consistent
          }else{               fixed_smk <- "tvSMK_"        }     # time-varying smoke
          
          
          mtype       <- ifelse(BaseModel,   "base_",    "primary_")
          if(SNPmainOnly){ mtype <- "SNPmain_" }
          fullresults <- ifelse(allresults,  "full",     "partial") 
          
                   
        #
          write.csv(all_gmmat, file=paste0(outfolder, "/", cohort, "_GMMAT_",  fixed_smk, mtype, fullresults, "_2024.csv"), row.names = F)   
          write.csv(all_gee,   file=paste0(outfolder, "/", cohort, "_GEE_",    fixed_smk, mtype, fullresults, "_2024.csv"), row.names = F)
          save(all_models,     file=paste0(outfolder, "/", cohort, "_summary_",fixed_smk, mtype, fullresults, "_2024.rdata") )
        
          write.csv(s_table,   file=paste0(outfolder, "/", cohort, "_table_2024.csv"), row.names = F)                           
        }
      
      

     return(all_models) 
     }







