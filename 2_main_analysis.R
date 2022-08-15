#----------Main analysis


#Overall model
rho <- 0.6
variance_matrix <- 
  impute_covariance_matrix(IB_dataset$vi, 
                           cluster = IB_dataset$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE) #compute varcov matrix

ml_model <- rma.mv(yi, V = variance_matrix, 
                   random = ~ 1 | StudyID/ESID,
                   data = IB_dataset) #main model
ml_model_robust <- robust(ml_model, cluster = IB_dataset$StudyID, clubSandwich = TRUE) #cluster robust adjustment
predict(ml_model_robust, transf=exp, digits = 2) #exponentiate to derive OR

#heterogeneity
i2 <- var.comp(ml_model) #calculation of I2; requires varcomp() function written by Mathias Harrer and David Ebert (Harrer, Cuijpers, Furukawa, & Ebert, 2019)
sum(ml_model_robust$sigma2) #tau^2


#-------Gender and Age
ml_model_age <- rma.mv(yi ~ age_imputed,
                       V = variance_matrix,
                       random = ~ 1 | StudyID/ESID ,
                       data=IB_dataset)
ml_model_age_robust <- robust(ml_model_age, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_age_robust, xlab="age", pred=TRUE, ci=TRUE)

ml_model_gender <- rma.mv(yi ~ male_prop_imputed,
                          V = variance_matrix,
                          random = ~ 1 | StudyID/ESID,
                          data=IB_dataset)
ml_model_gender_robust <- robust(ml_model_gender, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_gender_robust, xlab="proportion of males to females in sample", pred=TRUE, ci=TRUE)



#-------Literature
ml_model_lit <- rma.mv(yi ~ literature,
                       V = variance_matrix,
                       random = ~ 1 | StudyID/ESID,
                       data=IB_dataset)
ml_model_lit_robust <- robust(ml_model_lit, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_lit_robust, transf=exp, digits=2)



#---separate models for each literature (with dataset subsetted) are called for calculation of each theory's I2 and tau^2

#ml_model_AS <- rma.mv(yi, V = variance_matrix,
#                      random = ~ 1 | StudyID/ESID,
#                      data=IB_dataset,
#                      subset = (IB_dataset$literature == "AS"))
#ml_model_AS_robust <- robust(ml_model_AS, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
#i2 <- var.comp(ml_model_AS) calculation of I2; requires varcomp() function written by Mathias Harrer and David Ebert (Harrer, Cuijpers, Furukawa, & Ebert, 2019)
#sum(ml_model_AS_robust$sigma2)
#
#ml_model_load <- rma.mv(yi, V = variance_matrix,
#                        random = ~ 1 | StudyID/ESID,
#                        data=IB_dataset,
#                        subset = (IB_dataset$literature == "load"))
#ml_model_load_robust <- robust(ml_model_load, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
#i2 <- var.comp(ml_model_load) calculation of I2; requires varcomp() function written by Mathias Harrer and David Ebert (Harrer, Cuijpers, Furukawa, & Ebert, 2019)
#sum(ml_model_load_robust$sigma2)
