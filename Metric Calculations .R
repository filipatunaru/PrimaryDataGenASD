setwd("C:/Users/filip/OneDrive/Desktop/Roche Project/FILE COMPILE 2!!!!!!")

library(tidyverse)
library("ggplot2")
library("reshape2")
library(data.table)
library(readr)


memory.limit(9999999999)


data <- read_csv("data.csv")

data<- data # remove unnecessary columns at start

# SECTION 1 - DISTRIBUTION OF COVARIATES IN THE DATA SET

data.df <- as.data.frame(data)


dataQ72 <- data.df[data.df$Quarters_Since_Index == 72, ]

sum(dataQ72$EH_mental)/ nrow(data.df)
sum(dataQ72$EH_infection)/ nrow(data.df)
sum(dataQ72$EH_neoplasms)/ nrow(data.df)
sum(dataQ72$EH_metabolic)/ nrow(data.df)
sum(dataQ72$EH_blood)/ nrow(data.df)
sum(dataQ72$EH_neuro)/ nrow(data.df)
sum(dataQ72$EH_circ)/ nrow(data.df)
sum(dataQ72$EH_resp)/ nrow(data.df)
sum(dataQ72$EH_digestive)/ nrow(data.df)
sum(dataQ72$EH_genital)/ nrow(data.df)
sum(dataQ72$EH_pregnancy)/ nrow(data.df)
sum(dataQ72$EH_skin)/ nrow(data.df)
sum(dataQ72$EH_muscle)/ nrow(data.df)
sum(dataQ72$EH_congenital)/ nrow(data.df)
sum(dataQ72$EH_perinatal)/ nrow(data.df)
sum(dataQ72$EH_ill_def)/ nrow(data.df)
sum(dataQ72$EH_injury)/ nrow(data.df)



# proportion of regions of residence 

sum(dataQ72$east_midlands)/nrow(data.df)
sum(dataQ72$east_england)/nrow(data.df)
sum(dataQ72$london)/nrow(data.df)
sum(dataQ72$north_east)/nrow(data.df)
sum(dataQ72$north_west)/nrow(data.df)
sum(dataQ72$northern_ireland)/nrow(data.df)
sum(dataQ72$scotland)/nrow(data.df)
sum(dataQ72$south_central)/nrow(data.df)
sum(dataQ72$south_east_coast)/nrow(data.df)
sum(dataQ72$south_west)/nrow(data.df)
sum(dataQ72$wales)/nrow(data.df)
sum(dataQ72$west_midlands)/nrow(data.df)
sum(dataQ72$yorkshire_humber)/nrow(data.df)




# SECTION 2 - ASD DIAGNOSIS INFO PREVALENCE

# OVERALL PREVALENCE OF ASD

diag.per.quarter<- data.df %>% group_by(Quarters_Since_Index) %>% summarise(sum(ASD_diag))# cumulative

total.ASD <- max(diag.per.quarter[,2])
total.ASD
ASD.prev.data <- max(diag.per.quarter[,2])/ nrow(data.df)
ASD.prev.data

view(diag.per.quarter)




# MODAL AGE AT DIAGNOSIS 

newASD_diff <- diff(diag.per.quarter$`sum(ASD_diag)`)

newASD_diff <- c(0, newASD_diff) # number of new ASD diagnoses made per quarter 


diag.per.quarter$newASD_diff <- newASD_diff  


modal.quarters <- diag.per.quarter[diag.per.quarter$newASD_diff == max(diag.per.quarter$newASD_diff),]
modal.quarters


# MEDIAN AGE AT DIAGNOSIS 
median.info <- rep.int(diag.per.quarter$Quarters_Since_Index, diag.per.quarter$newASD_diff)
medianquarter <- median(median.info)   # expect this to be 24-27 (median age of diagnosis is 6) 
medianquarter


# DIAGNOSES BY AGE GROUP

#proportion of ASD diagnoses being made before age of 5 
ASD.diag.prop.05 <- diag.per.quarter[19,2] / total.ASD
ASD.diag.prop.05   # aim = 0.32

# proportion of ASD diagnoses being made between the ages of 5 and 12 
ASD.diag.prop.512 <- (diag.per.quarter[51,2]- diag.per.quarter[19,2])/ total.ASD
ASD.diag.prop.512 # aim = 0.58


# proportion of ASD diagnoses being made between age 13 and 18 
ASD.diag.prop.1318 <- (max(diag.per.quarter[72,2])- max(diag.per.quarter[51,2]))/ total.ASD

ASD.diag.prop.1318


# PLOT OF NEW DIAGNOSES OVER TIME
ggplot(diag.per.quarter,                  # Stacked barplot using ggplot2
       aes(x = Quarters_Since_Index ,
           y = newASD_diff)) +
  geom_col(stat = "identity")+ xlab( "Quarters Since Birth")+ 
  ylab("Number of New ASD Diagnoses")+ theme_classic()+theme(axis.title.x = element_text(
    colour="Black",size=10,face="bold"),
    
    axis.title.y = element_text(
      colour="Black",size=10,face="bold")
  )  + scale_y_continuous(expand = c(0, 0), limits = c(0, 430)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 72))  




#  SECTION 3 - SEX RATIOS IN SIMUALTED DATA 

# MALE TO FEMALE RATIO 
grouped.sex<- data.df %>% group_by(Patient_ID) %>% summarise(sum(sex))

grouped.sex

sum.sex <- sum(data.df$sex)/72

male.prop <- sum.sex / 600000
male.prop



# ASD MALE PROPORTION

df_ASD_Q72 <- data.df %>% filter(ASD_diag == 1 & Quarters_Since_Index == 72)

no.ASD.male <- sum(df_ASD_Q72$sex == 1, na.rm=TRUE) 

ASD.male.proportion <- no.ASD.male/ nrow(df_ASD_Q72)

ASD.male.proportion


# Distributions of new diagnoses over time 


data.female <- filter(data.df, sex == 0)

fem.diag<- data.female %>% group_by(Quarters_Since_Index) %>% summarise(sum(ASD_diag))# cumulative

new_fem_diag <- diff(fem.diag$`sum(ASD_diag)`)

new_fem_diag <- c(0, new_fem_diag)


diag.per.quarter$fem.new <- new_fem_diag

diag.per.quarter <- diag.per.quarter %>% mutate( male.new = newASD_diff - fem.new)


 


# plot stacked bar chart 

male.new.diag <- diag.per.quarter$male.new
fem.new.diag <- diag.per.quarter$fem.new


diag.data.sex <- data.frame(male.new.diag, fem.new.diag  )

sex.diag_long <- diag.data.sex              # Converting data from wide to long format
sex.diag_long$subgroup <- as.factor(1:nrow(sex.diag_long))
sex.diag_long <- melt(sex.diag_long, id.vars = "subgroup")
sex.diag_long

ggplot(sex.diag_long,                  # Stacked barplot using ggplot2
       aes(x = subgroup,
           y = value,
           fill = variable)) +
  geom_col(stat = "identity")+ xlab( "Quarters Since Birth")+ 
  ylab("Number of New ASD Diagnoses")+ theme_classic()+theme(axis.title.x = element_text(
    colour="Black",size=10,face="bold"),
    
    axis.title.y = element_text(
      colour="Black",size=10,face="bold")
  )  + scale_y_continuous(n.breaks = 15)+
  scale_x_discrete(breaks = seq(1, 72,2)) + labs (fill='Sex') 

