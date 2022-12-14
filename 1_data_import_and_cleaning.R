library(readxl)
library(dplyr)
library(ggplot2)
library(metafor)
library(clubSandwich)
library(PublicationBias)
library(robumeta)

#----------dataset import
IB_dataset <- read_excel("C:/Documents/IB_dataset.xlsx", na = "NA")

#----------data cleaning

#factors
factor_variables <- c('StudyID', 'literature', 'combined_theory', 'paradigm', 'setting', 'task', 'level', 'type', 'inhibition',
                      'load_type', 'PK', 'FA', 'TP', 'online', 'load_val', 'journal', 'country', 'o_l_manip', 
                      'age_range','task_imputed', 'q_score', 'published')

IB_dataset[,factor_variables] <- lapply(IB_dataset[,factor_variables], factor)

#quality score variable
IB_dataset$q_score <- ordered(IB_dataset$q_score, levels = c("100%_q", "66%_q", "33%_q", "0%_q"))

#gender proportion variable
IB_dataset$male_prop <- IB_dataset$males / (IB_dataset$males + IB_dataset$females)
IB_dataset$female_prop <- IB_dataset$females / (IB_dataset$males + IB_dataset$females)

#missing values for continuous variables
IB_dataset$male_prop_imputed <- replmiss(IB_dataset$male_prop, mean(IB_dataset$male_prop, na.rm = TRUE, digits=2))
IB_dataset$age_imputed <- replmiss(IB_dataset$age, mean(IB_dataset$age, na.rm = TRUE, digits=2))

#all duration missing data is from dynamic tasks
IB_dataset$CS_d <- replmiss(IB_dataset$CS_d, mean(IB_dataset$CS_d[IB_dataset$paradigm == "dynamic" & !IB_dataset$CS_d >= 15], na.rm = TRUE))
IB_dataset$trial_d <- replmiss(IB_dataset$trial_d, mean(IB_dataset$trial_d[IB_dataset$paradigm == "dynamic" & !IB_dataset$trial_d >= 60], na.rm = TRUE))
IB_dataset$CS_prop <- replmiss(IB_dataset$CS_prop, IB_dataset$CS_d / IB_dataset$trial_d)

#all target/distractor missing data is from attention set dynamic tasks; rounded because 3.5 targets/distractors cannot occur
IB_dataset$targets <- replmiss(IB_dataset$targets, round(mean(IB_dataset$targets[IB_dataset$literature == "AS" & IB_dataset$paradigm == "dynamic"], na.rm= TRUE)))
IB_dataset$distractors <- replmiss(IB_dataset$distractors, round(mean(IB_dataset$distractors[IB_dataset$literature == "AS" & IB_dataset$paradigm == "dynamic"], na.rm = TRUE)))


#------the following reverses the direction of effect sizes for cognitive load studies, ie reverse codes them, only unhash if running cog load reverse coded
#CL_studies <- IB_dataset$load_type == "cognitive"
#CL_studies[is.na(CL_studies)] = FALSE
#IB_dataset$yi[CL_studies] = -IB_dataset$yi[CL_studies]


#---outliers
hist(IB_dataset$yi)
mean_ES <- mean(IB_dataset$yi)
SD_ES <- sd(IB_dataset$yi)
outlier_bound <- SD_ES * 3.29
out_ub <- (mean_ES + outlier_bound)
out_lb <- (mean_ES - outlier_bound)
IB_dataset$outliers <- IB_dataset$yi >= out_ub | IB_dataset$yi<=out_lb

#a model is needed to calculate influential diagnostics (cooks distance, DF betas, hat values) >>>
rho <- 0.6
variance_matrix <- 
  impute_covariance_matrix(IB_dataset$vi, 
                           cluster = IB_dataset$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

ml_model <- rma.mv(yi, V = variance_matrix, 
                   random = ~ 1 | StudyID/ESID,
                   data = IB_dataset)


#---Influential diagnostics
#cooks distance
cd_IB <- cooks.distance.rma.mv(ml_model) #this takes a little while. grab a coffee!
plot(cd_IB, type="o", pch=19, ylab="Cook's Distance")
IB_dataset$cooks_distance <- cd_IB
IB_dataset$cooks_distance_outlier <- IB_dataset$cooks_distance > (4 / ml_model$k)
length(cd_IB)

sum(IB_dataset$cooks_distance_outlier)
which(IB_dataset$cooks_distance_outlier)
IB_dataset[IB_dataset$cooks_distance_outlier, c("Study", "ESID", "cooks_distance_outlier")]

#df betas
dfB_IB <- dfbetas(ml_model) #this also takes a little while. grab another coffee!
plot(dfB_IB$intrcpt, type="o", pch=19, ylab="DF betas")
IB_dataset$DFbetas <- dfB_IB
IB_dataset$DF_beta_outlier <- IB_dataset$DFbetas < -1 | IB_dataset$DFbetas > 1

sum(IB_dataset$DF_beta_outlier)
IB_dataset[IB_dataset$DF_beta_outlier, c("Study", "ESID", "DF_beta_outlier", "outliers")]

#hat values
hat_IB <- hatvalues(ml_model)
IB_dataset$hats <- hat_IB

k = as.numeric(ml_model$k)
p = as.numeric(length(ml_model$b))
hat_cutoff <- 3 * (p / k)
IB_dataset$hat_outliers <- IB_dataset$hats > hat_cutoff
sum(IB_dataset$hat_outliers)

# Note: descriptives are run prior to removing influentials/outliers

#remove outliers/influential cases from dataset
IB_dataset <- IB_dataset[!IB_dataset$outliers,]
IB_dataset <- IB_dataset[!IB_dataset$cooks_distance_outlier,]
