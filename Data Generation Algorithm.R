
library(tidyverse)
library(bindata)
library(LaplacesDemon)




task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)





##############################################################################
# DATA GENERATION
##############################################################################


total.no.patients <- 100000  # 100 batches of 100,000 patients were simulated
years.follow.up <- 18

no.quarters <- years.follow.up *4  # total number of quarters simulated

quarters.index <- seq(1, no.quarters, 1)  # quarters since index (birth)





# SIMULATING DISEASE COVARIATES 

no.covar <- 17 # number of risk factor variables, 17 disease categories which is also simulated from common 
#probability matrix. The time-invariant covariates of region and sex are simulated separately 


records <- matrix(ncol= no.covar  , nrow= no.quarters*total.no.patients) # create matrix for data to be stored in 


colnames(records)<- c("mental", "injury", "circ", "infection", "genital", "muscle", "pregnancy", "blood", "ill_def", 
                      "neuro", "neoplasms" , "metabolic", "digestive", "skin", "perinatal", "congenital", "resp")



covariate.prevalence <- c(0.048,0.281,0.011,0.503,0.135,0.112,0.005,0.017,0.304,
                          0.512,0.03,0.02,0.229,0.474,0.056, 0.063,0.648)  
# these have been calculated from the Roche case-control study and reflect the prevalence of the 17 disease categories in the population 
# across the 0-18 year follow up period 

covariate.prevalence <- covariate.prevalence/72  # above prevalence is for whole time period, since quarters are independently simulated 
# have to divide by 72 to maintain the above proportions in overall data 




# CREATING MATRIX OF COMMON PROBABILITIES FOR DISEASE COVARIATES

# The matrix of common probabilities:
mat <- outer(X = covariate.prevalence, Y = covariate.prevalence)

mat[lower.tri(mat)] <- NA #delete lower part, so that upper part can be transposed onto it later to make the matrix symmetric 

set.seed(task_id) # the seed should be set when generating full set of data 

mat <- mat * 2 ^ rnorm(n = length(mat), mean = 0, sd = 0.25) 

mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)] # to make matrix symmetric

diag(mat) <- covariate.prevalence  # replace diagonals with initial covariate prevalence

# effect sizes used for sex, disease categories and regions of residence 
betas <- c(1.1,0.5,0.2,0,-0.2,0,0,0,0,0.2,0,0,0.5,0.5,0.2,-0.2,0,0.2,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0.2)

# region prevalence
controls <- c(38437, 93083, 127035, 18104, 130716, 57878, 230006, 114813, 117653, 89141, 167296, 107948, 38276)

ASD_region <- c(185, 624, 737, 111, 860, 1556, 3657, 880, 1364, 495, 1701, 1002, 184)

regions.df <- data.frame(controls, ASD_region, row.names = c("east_midlands", "east_england", "london", "north_east", "north_west", "northern_ireland", 
                                                             "scotland","south_central", "south_east_coast", "south_west", "wales", "west_midlands", 
                                                             "yorkshire_humber"))

regions.df <- regions.df %>% mutate( region_total = ASD_region + controls, proportion = region_total/(1330386+13356), 
                                     no_control_notfromregion = 1330386- controls, no_ASD_notfromregion = 13356 - ASD_region)


region.prev <- regions.df$proportion


# baseline ASD risk function 
x <- seq(1,72,1)
beta.zero.func <- dlnorm(x, meanlog = log(25.5), sdlog = sqrt(log(25.5/17.5))) 
sum.female <- 0.2*sum(beta.zero.func)



# female scalars of beta.zero.func chosen such that baseline ASD diagnosis risk in females is concetrated in the later 
# quarters, but overall baseline female risk does not exceed 20% of overall baseline ASD risk

female.scalars <- c(  rep(0.0001, 10), 0.0002, 0.0003, 0.00015, 0.00016 ,0.00017, 0.00018,0.00019,0.0002, 0.00021,0.00022,# 1- 20
                      0.0025, 0.0027, 0.003, 0.0035, 0.004, 0.006, 0.009,0.01, 0.014, 0.015,#21-30
                      0.016, 0.018, 0.02, 0.022, 0.025, 0.027, 0.03, 0.035, 0.04 ,0.7, # 31-40 q. 44 = 12 years
                      0.8, 0.99, 0.99, 0.99, 0.99, 
                      0.99, 0.99, 0.99, 0.99, 0.99,  rep(0.99,5), rep(0.99,5 ), rep(0.99, 11) , 0)

female.diff <- sum.female - sum(beta.zero.func * female.scalars)

first.female.probs <- beta.zero.func * female.scalars

female.baseline <- c( first.female.probs[1:71], female.diff) # replace place holders with values 
# required so that sum of female baseline probabilities is 0.2*beta.zero.func (overall baseline probabilities)


male.baseline <- beta.zero.func - female.baseline 



baseline <- function(  patient.sex){
  if (patient.sex == 1) { 
    male.baseline
  }else {
    female.baseline
  } 
} 




###################################################################################################
# FUNCTION FOR GENERATING SINGLE PATIENT DATA 
#  - takes as arguments the number of quarters being generated, the matrix of common probabilities
#    and the regression coefficients 
###################################################################################################

patient.gen <- function(no.quarters, mat, beta ){
  
  covar <- suppressWarnings(rmvbin(n = no.quarters, commonprob = mat)) #creates correlated multivariate binary random variables 
  #by thresholding a normal distribution, 
  
  colnames(covar)<- c( "mental", "injury", "circ", "infection", "genital", "muscle", "pregnancy", "blood", "ill_def", 
                       "neuro", "neoplasms" , "metabolic", "digestive", "skin", "perinatal", "congenital", "resp")
  
  
  # GENERATING THE TIME-INVARIANT COVARIATES 
  
  # Sex 
  sex.df<- matrix(0, ncol =1, nrow = no.quarters)
  colnames(sex.df) = "sex"
  sex.draw <- rbinom(1,1, 0.51)
  sex.df <- as.data.frame(sex.df) %>% mutate(sex = case_when(sex.draw == 1 ~1, TRUE ~0))
  
  # Region
  reg.df <-  matrix(0,ncol= 13 , nrow= no.quarters) # create matrix for data to be stored in 
  colnames(reg.df) = c("east_midlands", "east_england", "london", "north_east", "north_west", "northern_ireland", 
                       "scotland","south_central", "south_east_coast", "south_west", "wales", "west_midlands", 
                       "yorkshire_humber")
  
  reg <- sample(1:13,1,replace=TRUE,prob=region.prev )# sample region from multinomial 
  
  reg.df <- as.data.frame(reg.df) %>% mutate(east_midlands = case_when(reg == 1 ~ 1, TRUE ~ 0),
                                             east_england = case_when(reg == 2 ~ 1, TRUE ~ 0),
                                             london = case_when(reg == 3 ~ 1, TRUE ~ 0),
                                             north_east = case_when(reg == 4 ~ 1, TRUE ~ 0),
                                             north_west = case_when(reg == 5 ~ 1, TRUE ~ 0),
                                             northern_ireland = case_when(reg == 6 ~ 1, TRUE ~ 0),
                                             scotland = case_when(reg == 7 ~ 1, TRUE ~ 0),
                                             south_central = case_when(reg == 8 ~ 1, TRUE ~ 0),
                                             south_east_coast = case_when(reg == 9 ~ 1, TRUE ~ 0),
                                             south_west = case_when(reg == 10 ~ 1, TRUE ~ 0),
                                             wales = case_when(reg == 11 ~ 1, TRUE ~ 0),
                                             west_midlands = case_when(reg == 12 ~ 1, TRUE ~ 0),
                                             yorkshire_humber = case_when(reg == 13 ~ 1, TRUE ~ 0))
  
  
  
  # Bind all covariate information together 
  covar2 <- cbind(sex.df, covar, reg.df )

  
  

  # ASD DIAGNOSIS GENERATION 
  
  eta <- as.matrix(covar2) %*% as.matrix( betas) %>% as.data.frame() %>% 
    mutate(Quarters_Since_Index = quarters.index, sex = sex.df)
  
  # select which baseline ASD risk values to use in calculation of ASD probability based on sex of generated patient
  baseline.values <- baseline( eta$sex[1,])
  
  
  ASD_diag <- rbinom(n = no.quarters, size = 1, prob = (invlogit(eta$V1-5.46+ 60*baseline.values[quarters.index]))/72)     

    single.patient <- cbind(covar2, ASD_diag) 
  
  for(i in 1:71){
    if(single.patient[i, "ASD_diag"] == 1){  #all entries for ASD_diag following first diagnosis (if occurs) are also changed to 1 
      single.patient[i+1, "ASD_diag"] <- 1
    }
  }
  

    # add quarter since index count
  single.patient <- single.patient%>% as.data.frame() %>%   mutate(Quarters_Since_Index = quarters.index, sex = sex.df)


  single.patient    # generates data for follow up period for a single person 

}


# Generate multiple patient records and collect in one dataframe 
patients <- 1:total.no.patients
dataout <- lapply(patients, 
                    function(s) patient.gen(no.quarters, mat, betas)) 
 
 all.patients.df <- do.call(rbind, dataout)

 
# Add patient ID 
  Patient.ID <- rep(c(seq(1, total.no.patients, 1)), each = no.quarters)
 all.patients.df <-  all.patients.df %>% mutate(Patient_ID = Patient.ID)

 
#########################################################################################
# MONITORING DATASET METRICS OF INTEREST
#########################################################################################
 
# Prevalence of ASD cases in all patients generated 

grouped.all <- all.patients.df %>% group_by(Patient_ID)%>% summarise(sum(ASD_diag))

num.ASD.data <- sum(grouped.all$`sum(ASD_diag)` >= 1, na.rm=TRUE)
ASD.prev.data <- num.ASD.data/ total.no.patients


# checking male:female ratio 

grouped.all.sex <- all.patients.df%>% group_by(Patient_ID) %>% summarise(sum(sex))
no.male <- sum(grouped.all.sex$`sum(sex)`>= 1, na.rm=TRUE)
male.proportion <- no.male/total.no.patients


# ASD male:female ratio

df_ASD_Q72 <- all.patients.df %>% filter(ASD_diag == 1 & Quarters_Since_Index == 72)
no.ASD.male <- sum(df_ASD_Q72$sex == 1, na.rm=TRUE)
ASD.male.proportion <- no.ASD.male/ nrow(df_ASD_Q72)


# checking distribution of ASD diagnoses over the time period 
diag.per.quarter<- all.patients.df %>% group_by(Quarters_Since_Index) %>% summarise(sum(ASD_diag))


# proportion of ASD diagnoses being made before age of 5 
ASD.diag.prop.05 <- (max(diag.per.quarter[1:19,2])) / num.ASD.data
 aim = 0.32

# proportion of ASD diagnoses being made between the ages of 5 and 12 
ASD.diag.prop.512 <- (max(diag.per.quarter[51,2])- max(diag.per.quarter[19,2]))/ num.ASD.data
 # aim = 0.58


# proportion of ASD diagnoses being made between age 13 and 18 
ASD.diag.prop.1318 <- (max(diag.per.quarter[72,2])- max(diag.per.quarter[51,2]))/ num.ASD.data
  # aim = ~ 0.09


# number of new ASD diagnoses per quarter 
newASD_diff <- diff(diag.per.quarter$`sum(ASD_diag)`)

newASD_diff <- c(0, newASD_diff)

diag.per.quarter$newASD_diff <- newASD_diff


# Modal age of ASD diagnosis
modal.quarters <- diag.per.quarter[diag.per.quarter$newASD_diff == max(diag.per.quarter$newASD_diff),]
modalquarter <- modal.quarters$Quarters_Since_Index# would expect this to be 16-19 (modal age of diagnosis is 4) 


# median age of ASD diagnosis
median.info <- rep.int(diag.per.quarter$Quarters_Since_Index, diag.per.quarter$newASD_diff)
medianquarter <- median(median.info)   # expect this to be 24-27 (median age of diagnosis is 6) 




##############################################################################################
# Creating "Summary Check" Table to Display Metrics Calculated Above
##############################################################################################

check.names <- c("ASD_prev", "Total_Male", "ASD_Male", "DIAG_<5", "DIAG_5-12", "DIAG_13_18",
                  "Median.Quarter")

results <- c(ASD.prev.data, male.proportion, ASD.male.proportion, ASD.diag.prop.05, ASD.diag.prop.512, 
            ASD.diag.prop.1318,  medianquarter )

summary.checks <- data.frame(check.names, results)


# Add column to check whether data metric are within 10% of target 
summary.checks[1, 3] <- case_when( 0.0165 < summary.checks[1, 2] &  summary.checks[1, 2] <= 0.021 ~ "YES", TRUE ~ "NO")
summary.checks[2, 3] <- case_when( 0.4845 < summary.checks[2, 2] &  summary.checks[2, 2] <= 0.5355 ~ "YES", TRUE ~ "NO")
summary.checks[3, 3] <- case_when( 0.7505 < summary.checks[3, 2] &  summary.checks[3, 2] <= 0.8295 ~ "YES", TRUE ~ "NO")
summary.checks[4, 3] <- case_when( 0.288 < summary.checks[4, 2] &  summary.checks[4, 2] <= 0.352 ~ "YES", TRUE ~ "NO")
summary.checks[5, 3] <- case_when( 0.522 < summary.checks[5, 2] &  summary.checks[5, 2] <= 0.638 ~ "YES", TRUE ~ "NO")
summary.checks[6, 3] <- case_when( 0.081 < summary.checks[6, 2] &  summary.checks[6, 2] <= 0.099 ~ "YES", TRUE ~ "NO")
summary.checks[7, 3] <- case_when( 24 < summary.checks[8, 2] &  summary.checks[8, 2] <= 27 ~ "YES", TRUE ~ "NO")




#############################################################################################
# APPENDING AUGMENTED DISEASE CATEGORY VARIABLES TO DATA
#############################################################################################

# Appending Cumulative Disease Category Variables

df_aug <- all.patients.df %>% group_by(Patient_ID)%>% mutate(cum_mental = cumsum(mental),
                                      cum_injury = cumsum(injury),
                                      cum_circ = cumsum(circ),
                                      cum_infection = cumsum(infection),
                                      cum_genital = cumsum(genital),
                                      cum_muscle = cumsum(muscle),
                                      cum_pregnancy = cumsum(pregnancy),
                                      cum_blood = cumsum(blood),
                                      cum_ill_def = cumsum(ill_def),
                                      cum_neuro= cumsum(neuro),
                                      cum_neoplasms = cumsum(neoplasms), 
                                      cum_metabolic = cumsum(metabolic),
                                      cum_digestive = cumsum(digestive),
                                      cum_skin = cumsum(skin),
                                      cum_perinatial = cumsum(perinatal),
                                      cum_congenital = cumsum(congenital),
                                      cum_resp = cumsum(resp)) 


# Appending "Ever-Had" Disease Category Variables 

df_EH <- df_aug[,35:51] # subset the cumulative disease columns, so they can be transformed into "ever-had" columns 

df_EH<-  replace(df_EH, df_EH >=1, 1)
colnames(df_EH) <- c( "EH_mental", "EH_injury", "EH_circ", "EH_infection", "EH_genital", "EH_muscle", "EH_pregnancy", 
                      "EH_blood", "EH_ill_def", "EH_neuro", "EH_neoplasms" , "EH_metabolic", "EH_digestive", "EH_skin", 
                      "EH_perinatal", "EH_congenital", "EH_resp")

data <- cbind(df_aug, df_EH)  

data <- data  %>% select(Patient_ID, sex, Quarters_Since_Index, everything()) # reorder columns for ease of reading data



# Saving Data and Summary Check 

filename.data <- paste0("datafile_", task_id, ".RData")

filename.summary <- paste0("summaryfile_", task_id, ".RData")

 
save(data , file= filename.data) 


save(summary.checks, file= filename.summary)





