
###############################################################
#   use (age-20) to simulate fev1 & add steps for HWE
#    generate "pure" slope and baseline affect for SNP
###############################################################
# module load  R/4.0.2


##################################################################
#  NOTE ON INPUT DATA
# simulation data used in manuscript was generated based on the real FHS age, which is not publicly available.



## =====================================================
## effect size (e3)  
   library(MASS)
  
## common parameters
   N          <- length(unique(d$id))
   decline_pc <- 0.1    # decline will be 10% of the main effect of AGE-20 
   snp_pc     <- 0.1    # interaction effect of SNP*(age-20)

## covariance
   Y_UN1_1    <- 0.18
   Y_UN2_1    <- -0.00064 
   Y_UN2_2    <- 0.000131
   error_var  <- 0.021
   AR1        <- 0.7
   mcov       <- matrix(c(Y_UN1_1, Y_UN2_1, Y_UN2_1, Y_UN2_2), nrow = 2)
   SD_resid        <- sqrt(error_var)
   SD_resid_unique <- sqrt(1-AR1*AR1)*SD_resid
  


# ========================================================
  f_simdata <- function(d, Ngroup=7, LorNL, snp_effect, nullmodel=FALSE){ # Ngroup=4; LorNL=1; snp_effect=1
      # coefficients
        if(LorNL == 1){
            mean_intercept <- 3
            mean_b1        <- -0.018
            
        }else if(LorNL == 2){
            mean_intercept <- 3.25
            mean_b1        <- -0.0115
            fixed_b2       <- -0.00019           
        }else{
            print("No valid input")
        }
        mu <- c(mean_intercept, mean_b1)
      
     # -------------------------------------------------------
       ids  <- unique(d$id)
       bi   <- mvrnorm(n = N, mu, mcov, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
       if(LorNL == 1){      
          dat           <- cbind(ids, bi)                               #linear   
          colnames(dat) <- c("ID", "intercept", "B1")          
       }else if(LorNL == 2){  
          dat           <- cbind(ids, bi, rep(fixed_b2, length(ids)))   #nonlinear   
          colnames(dat) <- c("ID", "intercept", "B1", "B2")          
       }else{print("Invalid number")}
      
       dat <- merge(d, dat, by.x="id", by.y="ID", all.x = T) 
                
     # ------------------------------------------------------- 
     # generate residuals
       dat <- dat[order(dat$id, dat$timefactor_spiro), ]   # order matters !!!!
       dat$residual <- rep(0, nrow(dat))
       for(j in 1:nrow(dat)){  
           if(dat$timefactor_spiro[j] == 0){ 
              dat$residual[j] <- rnorm(1, mean=0, sd=SD_resid)
           }else{ 
              dat$residual[j] <- AR1^(dat$age[j]-dat$age[j-1])*dat$residual[j-1] + rnorm(1, mean=0, sd=SD_resid_unique) 
           }
         }
    
     # -------------------------------------------------------- 
     # generate outcome (FEV1) with NORMAL decline
     # Linear
       if(LorNL == 1){
             dat$fev1 <- dat$intercept + dat$B1*dat$age_20 + dat$residual 
     # nonLinear   
       }else if(LorNL == 2){   
             dat$fev1 <- dat$intercept + dat$B1*dat$age_20 + dat$B2*(dat$age_20^2) + dat$residual 
       }else{print("Invalid number")}
  
     
     # add number of PFT count
       n_pft           <- data.frame(table(dat$id))
       colnames(n_pft) <- c("id", "n_pft")
       dat             <- merge(dat, n_pft, by = "id", all.x = T)
    
  
     ## =================================================================
     ## generating groups based on SNP effects using baseline data
     ##   snp_effect 
     ##   1 = snp affects decline(slope) 
     ##   2 = snp affects baseline
     ##   3 = no snp effects
     ## =================================================================     
     
       dat_base <- dat[which(dat$timefactor_spiro==0), c("id", "n_pft", "age")] # 3258
       colnames(dat_base)[3] <- "age_base"
     
     # Effect allele frequency 
       p      <- 0.2             # Since we want to study common variant, MAF should > 0.05
       freq_0 <- (1 - p)^2       
       freq_1 <- 2 * p * (1 - p) 
       freq_2 <- p^2             

      ## ================================================================
      ## 1.  Decline-type classes present in the population
      ##     Sample decline type  (independent of genotype) 
      ## ================================================================
         groups      <- c("normal", "lb1", "lb2", "rd1", "rd2", "lbrd1",  "lbrd2")
         group_prob  <- c( normal = 0.6,
                           lb1    = 0.1,   lb2    = 0.05,
                           rd1    = 0.1,   rd2    = 0.05,
                           lbrd1  = 0.05,  lbrd2  = 0.05 )
      
         dat_base$grp <- sample(groups, N, replace = TRUE, prob = group_prob)
      
      ## ================================================================
      ## 2.  Sample SNP genotype under HWE  (independent of grp)
      ## ================================================================
         dat_base$snp_true <- sample(c(0, 1, 2), N, replace = TRUE, prob = c(freq_0, freq_1, freq_2))
         dat_base$snp      <- dat_base$snp_true
      
      ## ================================================================
      ## 3. impose the *causal* SNP effect
      ##    If nullmodel == TRUE, skip this section.
      ## ================================================================
         b_rd          <- mean_b1 * decline_pc       
         b_lowbaseline <- 0.25                       
         snp_b_int     <- (mean_b1+b_rd)/2 * snp_pc  
         snp_base      <- -0.025                     

         if(!nullmodel){
           if(snp_effect == 1){         dat_base$b_int      <- snp_b_int   
           }else if(snp_effect == 2){   dat_base$delta_base <- snp_base     
           }
         }

      ## ================================================================
      ## ================================================================
         dat <- merge(dat, dat_base, by = c("id", "n_pft"))
        
      ## baseline shifts by grp 
         dat$fev1 <- ifelse(dat$grp %in% c("lb1", "lbrd1"), dat$fev1 - b_lowbaseline/2, dat$fev1)
         dat$fev1 <- ifelse(dat$grp %in% c("lb2", "lbrd2"), dat$fev1 - b_lowbaseline,   dat$fev1)
        
      ## slope shifts by grp 
         dat$fev1 <- ifelse(dat$grp %in% c("rd1", "lbrd1"), dat$fev1 + (b_rd/2) * dat$age_20, dat$fev1)
         dat$fev1 <- ifelse(dat$grp %in% c("rd2", "lbrd2"), dat$fev1 + b_rd * dat$age_20,     dat$fev1)
        
      ## add the SNP effect *after* baseline/slope tweaks
         if(!nullmodel){
           if(snp_effect == 1){
               dat$fev1 <- dat$fev1 + dat$snp_true * snp_b_int * dat$age_20
           }else if(snp_effect == 2){
               dat$fev1 <- dat$fev1 + dat$snp_true * snp_base
           }
         }
      
     # add additional variables
       dat$time         <-  dat$age - dat$age_base
       dat$timesq       <-  dat$time^2
       dat$baseage_time <-  dat$age_base * dat$time
       ##############################################
       dat$age          <- (dat$age - 20)
       dat$ageC2        <- (dat$age - 20)^2
       dat$ageXtime     <-  dat$age * dat$time
       ############################################## 
     # Replace negative FEV1 values by 0
       dat$fev1 <- ifelse(dat$fev1<0, 0, dat$fev1) 
  
  
     final <- list("dat"=dat)
     return(final)
  }
