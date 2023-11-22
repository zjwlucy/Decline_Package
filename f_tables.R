## functions to 
## generate summary tables 


## ---------------------------------------------------------------
## Table 1: Characteristics of study participants
# Columns: characteristic, study name 1, study name 2…..

# Rows:  n (observations), age, sex, race, height, pack years, smoking status , FEV1/FVC, FEV1 % predicted, FEV1, number of observations (total [ min, max]).  


# For continuous variable:   mean (sd)
# For categorical variable:  count of variable = 1  (%)

     
# baseline characteristics of data     
  f_tables <- function(dat, multiRace){
    print("----------Generating summary table for the full dataset----------")
    print("Note: variables must be coded in the format shown in instruction")
       
    tryCatch({      
         
      ## variables need to be summarized 
         dat      <- dat[order(dat$FID, dat$IID, dat$timefactor_spiro),]
         var_cont <- c("pre_fev1", "pre_fev1fvc", "fev1_pp", "age", "ht_baseline", "smoking_packyears_base")
         var_cat  <- c("sex", "smoking_status_base")  # "smoking_status"
         
         if(multiRace){   var_cat <- c(var_cat, "race")   }
         
      ## Take variables that are available from the data given    
         var_cont <- base::intersect(var_cont, colnames(dat))
         var_cat  <- base::intersect(var_cat,  colnames(dat))
         
         
      ## add Rank for observations 
         dat_count <- lapply(unique(dat$IID), function(x){
                             tmp         <- dat[which(dat$IID == x), c("IID", "obsID", "age")]
                             tmp$obsRank <- rank(tmp$age)
                             return(tmp)
                             })
         dat_count <- do.call(rbind, dat_count)
         dat       <- merge(dat, dat_count[,c("IID", "obsID", "obsRank")], by=c("IID", "obsID"), all.x=T)                   
                 
      ## add n_pft
         pft_count           <- as.data.frame(table(dat$IID))    # add number of pfts for each individual
         colnames(pft_count) <- c("IID", "n_pft")                # rep(pft_count$n_pft, pft_count$n_pft)
         dat                 <- merge(dat, pft_count, by="IID", all.x=T)
      
         
            
      ##----------------------------------------------------------        
      # Sample size:  total obs (total unique individuals)
        a_n <- paste0(nrow(dat), " (", length(unique(dat$IID)), ")")
        
      # Time: average of the total follow up duration: (Last-Baseline)=max(timefactor_spiro)  # d2 <- dat %>%group_by(IID)%>%slice(n()) 
        d_follow <- lapply(unique(dat$IID), function(x){
                           tmp  <- max(dat$timefactor_spiro[which(dat$IID == x)])
                           tmax <- c(x, tmp)
                          })   # 3258
        d_follow <- as.data.frame(do.call(rbind, d_follow))
        colnames(d_follow) <- c("IID", "t")
         
        a_time <- rbind(
                  # follow up time for ALL individuals
                    paste0( round(mean(d_follow$t, na.rm = T), digits = 1), 
                        " (", round(sd(d_follow$t,na.rm = T), digits = 1), ") "),
                  # follow up time excluding individuals with only 1 visit      
                    paste0( round(mean(d_follow$t[which(d_follow$t>0)], na.rm = T), digits = 1), 
                        " (", round(sd(d_follow$t[which(d_follow$t>0)],na.rm = T), digits = 1), ") ")  
                        )
        
        a_n    <- as.data.frame(a_n)
        a_time <- as.data.frame(a_time)
       
        a_n$variable    <- "N"
        a_time$variable <- c("TotalFollowUp", "TotalFollowUp_exclude_IndWith1Visit") 
                
        colnames(a_n)    <- c("value", "variable")
        colnames(a_time) <- c("value", "variable") 
       
        a_n$types    <- "sample_size_Total(Unique)"
        a_time$types <- rep("continuous_BaseToLast",2)
       
       
      ##----------------------------------------------------------- 
      # baseline characteristics 
        dat_base  <- dat[which(dat$obsRank == 1), ]

        #
          a_pft <- paste0(round(mean(dat_base$n_pft,na.rm = T), digits=0) , 
                           " [", min(dat_base$n_pft,na.rm = T),  
                           ",",  max(dat_base$n_pft,na.rm = T), "]")
          a_pft           <- as.data.frame(a_pft)
          a_pft$variable  <- "N_pft"  
          colnames(a_pft) <- c("value", "variable")
          a_pft$types     <- "count_pft(min,max)"              
           
        #                                                         
          a_cont <- apply(dat_base[, var_cont], 2, function(x)  
                          paste0(round(mean(x,na.rm = T), digits = 1), " (", round(sd(x,na.rm = T), digits = 1), ")") )                                                   
          a_cont           <- as.data.frame(a_cont)
          a_cont$variable  <- rownames(a_cont)
          colnames(a_cont) <- c("value", "variable") 
          a_cont$types     <- rep("continuous_base", nrow(a_cont))
          
        # pack years among smokers
          a_pkyrs <- apply(as.data.frame(dat_base[which(dat_base$smoking_status_base!=0), c("smoking_packyears_base")]), 
                           2, function(x)  
                           paste0(round(mean(x,na.rm = T), digits = 1), " (", round(sd(x,na.rm = T), digits = 1), ")") )
          names(a_pkyrs) <- "smoking_packyears_base_smokers"
          a_pkyrs        <- data.frame(value=a_pkyrs, variable="smoking_packyears_base(among smokers)", types="continuous_base")
          a_cont         <- rbind(a_cont, a_pkyrs)

            
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
