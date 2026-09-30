################################################
#       creating example simulated data        #                     
################################################
# module load  R/4.0.2  
  

############################ Note ##################################
# Simulation requires "your_longitudinal_data.csv"
# which is FHS real data and not publicly avaiable. 
# We provided a "fake_data_n1000.csv" with 100 randomly generated individuals as an example, for test run

#  "your_longitudinal_data.csv": long-format data containing individuals as rows and has an "age" column,
# where each individual may have multiple rows.
#  To run this script, supply your own long-format data
#  with these columns:
#    id               : participant ID
#    exam             : exam number (used to order visits)
#    count            : visit number within participant (1, 2, ...)
#    timefactor_spiro : 0 at the baseline visit, > 0 at follow-ups
#    age              : age (years) at each visit
#  Each participant must have exactly one baseline row (timefactor_spiro == 0)
##################################################################



# For 7 groups, snp effects slope
  dset        <- 2    # 1 for linear decline; 2 for nonlinear decline
  snp_effect  <- 1    # 1 = snp affects decline(slope);  2 = snp affects baseline; 3 = no effect   

  
# Data for age distribution
  #d <- read.csv("your_longitudinal_data.csv", na.strings="")  
  d           <- read.csv("fake_data_n1000.csv", na.strings="")  
  d           <- d[order(d$id, d$exam), ]   
  rownames(d) <- NULL
  d$age_20    <- d$age - 20
  
  
# 
  source("generate_data_v5.R")

  set.seed(2026)
  d_all   <- f_simdata(d=d, LorNL=dset, snp_effect=snp_effect)
  dat_all <- d_all$dat
  dat_all <- dat_all[order(dat_all$id, dat_all$count, dat_all$age),]
  rownames(dat_all) <- NULL
 
  dset <- ifelse(dset == 1, "L", "nonL")
  dat_all$dset       <- dset
  dat_all$Ngroup     <- 7
  dat_all$snp_effect <- snp_effect
  
  
  save(dat_all, file=paste0("d_", dset, "_", snp_effect, ".Rdata"))
         
         
 
