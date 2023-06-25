## functions to 
## generate summary tables 


## ---------------------------------------------------------------
## Table 1: Characteristics of study participants
# Columns: characteristic, study name 1, study name 2�..

# Rows:  n (observations), age, sex, race, height, pack years, smoking status , FEV1/FVC, FEV1 % predicted, FEV1, number of observations (total [ min, max]).  


# For continuous variable:   mean (sd)
# For categorical variable:  count of variable = 1  (%)

     
# baseline characteristics of data     
  f_tables <- function(dat, multiRace){
    print("----------Generating summary table for the full dataset----------")
    print("Note: variables must be coded in the format shown in instruction")
       
    tryCatch({      
         
      ## variables need to be summarized 
         dat      <- dat[order(dat$IID, dat$timefactor_spiro),]
         var_cont <- c("pre_fev1", "pre_fev1fvc", "fev1_pp", "age", "smoking_packyears_base", "ht_baseline")
         var_cat  <- c("sex", "smoking_status")
       
         if(multiRace){
            var_cat <- c(var_cat, "race")
         }
            
            
      ##----------------------------------------------------------        
      # sample size:  total obs (total unique individuals)
      # time: needs to be summarized using the longitudinal data not the baseline data       
        a_n    <- paste0(nrow(dat), " (", length(unique(dat$IID)), ")")
        a_time <- paste0( round(mean(dat$timefactor_spiro, na.rm = T), digits = 1), " (", round(sd(dat$timefactor_spiro,na.rm = T), digits = 1), ") ") 
       
        a_n    <- as.data.frame(a_n)
        a_time <- as.data.frame(a_time)
       
        a_n$variable    <- "N"
        a_time$variable <- "time"
        
        colnames(a_n)    <- c("value", "variable")
        colnames(a_time) <- c("value", "variable") 
       
        a_n$types    <- "sample_size"
        a_time$types <- "continuous_long" 
       
      ##----------------------------------------------------------- 
      # baseline characteristics  
        dat_base <- dat[which(dat$timefactor_spiro == 0), ]
         
      # check if number of observations match
        pft_count            <- as.data.frame(table(dat$IID))      # add number of pfts for each individual
        colnames(pft_count)  <- c("IID", "n_pft")   
        if( identical(dat$n_pft, rep(pft_count$n_pft, pft_count$n_pft)) ){
        
          #
            a_pft <- paste0(round(mean(dat_base$n_pft,na.rm = T), digits=0) , 
                            " [", min(dat_base$n_pft,na.rm = T),  
                            ",",  max(dat_base$n_pft,na.rm = T), "]")
            a_pft           <- as.data.frame(a_pft)
            a_pft$variable  <- "N_pft"  
            colnames(a_pft) <- c("value", "variable")
            a_pft$types     <- "count"              
           
          #                                                         
            a_cont <- apply(dat_base[, var_cont], 2, function(x)  
                                paste0(round(mean(x,na.rm = T), digits = 1), " (", round(sd(x,na.rm = T), digits = 1), ") ") )                                                              
            a_cont <- as.data.frame(a_cont)
            a_cont$variable  <- rownames(a_cont)
            colnames(a_cont) <- c("value", "variable") 
            a_cont$types     <- rep("continuous_base", nrow(a_cont))
            
          #    
            a_cat  <- apply(dat_base[, var_cat], 2, function(x) {
                                    test      <- as.data.frame(table(x))
                                    test$perc <- round(test$Freq/sum(test$Freq), digits=2) 
                                    return(test)} )              
            a_cat  <- do.call(rbind, a_cat)
            a_cat$variable  <- rownames(a_cat)
            a_cat$variable  <- gsub("\\..*", "", a_cat$variable)
            a_cat$variable  <- paste0(a_cat$variable, "_", a_cat$x)
            a_cat$value     <- paste0(a_cat$Freq, " (", a_cat$perc*100, "%)" )
            rownames(a_cat) <- NULL
            a_cat           <- a_cat[, c("value", "variable")] 
            a_cat$types     <- rep("categorical_base", nrow(a_cat))
          
            s_all <- rbind(a_n, a_pft, a_time, a_cont, a_cat)
            rownames(s_all) <- NULL
            s_all <- s_all[, c("variable", "value", "types")]
            
               
         }else{ print("ERROR: missmatch in the number of observations for each individual")}
           

       }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
       
  return(s_all)
  }






## ------------------------------------------------------------
## Table 2: Simulations - pick top LMM and GEE models and compare to reference model (two time point model). Remainder of models in supplement

# Columns, model, interaction term(s), random effects, null SNP (avg beta, avg p, % best model, % converged), low baseline SNP (avg beta, avg p, % best model, % converged), rapid decline SNP (avg beta, avg p, % best model, % converged).
# Rows: age + time models 







## ------------------------------------------------------------
## Table 3: Top clinical covariate model in real cohorts + snp effects; compare to reference model of 2 time point model

# Rank models by clinical covariate model. Show snp X time variable results for each cohort. Calculate avg AIC, avg computational time, and avg % converged.
# Columns: models, covariates, FHS (snp X time, p), EC (snp X time, p), COPDGene (snp X time, p),...pooled cohorts (snp X time, p), avg AIC, avg computational time, avg % converged.
# Rows: use covariate model as first row showing effect size of null snp X time, then show individual SNPs added to these.