# Decline_Package

## 1. Background
The optimal approaches to study lung function decline, particularly in the setting of omic-scale predictors, are not known. The goal of this research project is to test different statistical models (e.g. slope, GEE, GLM) in simulated and real world cohorts (population and case-control) to compare Type 1 and Type 2 error, identify factors leading to heterogeneity, and to provide a foundation for modeling for future omics data.

### Aims
1. develop a population of simulated spirometry data with a) autocorrelation, estimated from real world data, b) four groups representing normal, lower baseline, rapid decline, and lower baseline and rapid decline and c) both linear and quadratic decline
2. (this code package): test a set of models in real world data: a) simple slope models; b) LMM with random intercept, slope, or both; c) GEE; d) age, time; e) quadratic or linear terms (for underlying linear or quadratic decline). We will look at the effects of smoking (as a positive control), and selected SNPs.

### Authorship
Anticipate up to 4 per cohort, with additional as required for writing / analysis.

### Methods
Prepare cohort dataset as requested below. We require pre_fev1 in mL, age, sex, height, race, smoking status, pack-years,  SNP allele freq, fev1pp (can use GLI global, or can use what has been previously calculated for your cohort), fev1/fvc ratio;  follow up time, number of visits. We anticipate that the dataset will be clean; i.e. with minimal missingness, subjects with existing longitudinal data (>= 2 time points, smoking data), removal of erroneous data (i.e. QC’d spirometry and identification of spurious outliers). Large discrepancies in sample size in models / baseline characteristics will be assessed, and if necessary, request to re-prepare the datasets. We selected 12 SNPs: 6 from prior GWAS (COPD, lung function, or lung function decline) and 6 proxies (null control). 
If you have related individuals, please reach out to Jingwen.
Errors particularly with GEE model 13 are known, the software will continue to run. 

### Results
All result files will be added to the source folder.
Please zip/tar/etc with filename cohort_date.zip and send to Jingwen (zjwlucy@bu.edu) and Matt Moll (remol@channing.harvard.edu)
See project document for further details on the project.


## 2. Columns for data set 
  * __`IID`__:                    Unique individual ID (numeric). Note: two different individuals CANNOT have the same IID
  * __`FID`__:                    family ID (numeric); For unrelated individuals, create a column of the same number as "FID", e.g. a column of 1
  * __`pre_fev1`__:               FEV1 (in mL)          
  * __`SNPs`__:                   SNP information, column name MUST starts with the lower case "rs", e.g. rs507211, rsChrPosRefAlt
  * __`timefactor_spiro`__:       time since the baseline exam (in YEARS). At baseline, timefactor_spiro=0
  * __`age`__:                    time-varying age
  * __`smoking_status`__:         time-varying smoking status (never=0, former=1, current=2); 

### Baseline 
  * __`age_baseline`__:           baseline age
  * __`ht_baseline`__:            baseline height (in cm)
  * __`smoking_packyears_base`__: pack-years at baseline
  * __`sex`__:                    biological sex (female=0, male=1)
  * __`smoking_status_base`__:    baseline smoking status (never=0, former=1, current=2); Will be used as the grouping variable for glmmkin.
     
### Other cohort-specific variables
  * PCs (PC1, PC2, ... ...) 
  * equipchange ......  

### Variables for summary not for analysis
  * __`pre_fev1fvc`__:   ratio of fev1 and fvc (fev1/fvc) as percentage (e.g., fev1/fvc=0.75, the value for this column will be 75)
  * __`fev1_pp`__:       fev1 percent predicted (e.g., 92% as 92 not 0.92) 


###  
  * Kinship matrix (for related data):   both row names and column names MUST be IID      
      
## 3. Important Note
  * Cohorts with multiple racial groups should conduct **race-stratified** analyses.
  * If cohorts have only 3 or less repeated measurements per individual, you may encounter an error message for **GEE model 13** (e.g. something related to contrasting for variables with less than 2 levels). If this error only occurs in GEE model 13, please ignore this error, but kindly notify us and provide results for the remaining models.
  * Please check the files from the output folder "decline_package_output" to make sure that you have all the requested files and plots (plots are displayed normally). 

## 4. Required packages
If the current version does not work, you can try previous version.
Both versions should give the SAME output

### R version: R/4.2.1   (current version)
   * geepack_1.3.9
   * GMMAT_1.4.0 
   * dplyr_1.1.2  
   * readxl_1.4.3
   * ggplot2_3.4.2     
   * kinship2_1.9.6 (optional) for kinship matrix 

### R version: R/4.0.2  (previous version)
   * geepack_1.3-2
   * GMMAT_1.3.1
   * dplyr_1.0.2
   * readxl_1.3.1
   * ggplot2_3.3.2
   * kinship2 (optional) 

## 5. For analysis, use "Analysis_MAIN.R"
Example code is included inside Analysis_MAIN.R as comments. Please make sure your dataset has all the columns that are mentioned in **section 2. Columns for data set**.  

Example data: 
|IID	|FID|	pre_fev1|	timefactor_spiro|	age|	smoking_status|age_baseline|	…	|PC1|rs507211|	rs4077833|
|----|---|---------|-----------------|----|---------------|------------|---|-----|--------|----------|
|1	|1	 |2771	|0	 |50	|1	|50	|…	|-0.0293|0	|1|
|1	|1	 |2500	|4	 |54	|1	|50	|…	|-0.0293|0	|1|
|1	|1	 |2450	|10	|60	|1	|50	|…	|-0.0293|0	|1|
|2	|11	|3510	|0	 |38	|0	|38	|…	|-0.0038|0	|2|
|2	|11	|3450	|9	 |47	|0	|38	|…	|-0.0038|0	|2|
|2	|11	|3320	|12	|50	|0	|38	|…	|-0.0038|0	|2|
|2	|11	|3220	|17	|55	|0	|38	|…	|-0.0038|0	|2|
|3	|24	|2570	|0	 |42	|2	|42	|…	| 0.0071|1	|0|
|3	|24	|2600	|6	 |48	|2	|42	|…	| 0.0071|1	|0|
|3	|24	|2540	|8	 |50	|2	|42	|…	| 0.0071|1	|0|


Example Kinship matrix for 5 individuals (IID = 11,20,31,42,50):
|IID |11	|20	|31	|42	|50|
|---|-----|-----|-----|-----|-----|
|**11**	|0.5	  |0	    |0.25	 |0.125	|0     |
|**20**	|0	    |0.5	  |0.25	 |0.125	|0     |
|**31**	|0.25 	|0.25	 |0.5	  |0.25	 |0.25  |
|**42**	|0.125	|0.125	|0.25	 |0.5	  |0     |
|**50**	|0	    |0	    |0.25	 |0	    |0.5   |

