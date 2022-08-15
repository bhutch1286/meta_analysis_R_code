#----------Inhibition analysis


#1.Inhibition studies versus not inhibition studies
inhib_v_no_inhib <- rma.mv(yi ~ inhibition,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset,
                           subset = (IB_dataset$literature == "AS"))
inhib_v_no_inhib_robust <- robust(inhib_v_no_inhib, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
predict(inhib_v_no_inhib_robust, transf=exp, digits=2)


#2. Match versus neutral and match versus mismatch
inhibition_analysis <- read_excel("C:/Documents/inhibition.xlsx")

#remove neutral vs inhibition comparisons
inhibition_analysis$matchv1 <- inhibition_analysis$comparison != "neutral_v_inhibition"
inhibition_analysis_match_test <- inhibition_analysis[inhibition_analysis$matchv1,]
#inhibition_analysis_match_test <- inhibition_analysis[inhibition_analysis$matchv1 & !inhibition_analysis$low_prec_removed,] #unhash for low precision study trim analysis

#low precision cases are removed studywise in the inhibition analysis because the inhibition analysis is made up of different comparisons than the main analysis

rho <- 0.6
variance_matrix_inhibition <- 
  impute_covariance_matrix(inhibition_analysis_match_test$vi, 
                           cluster = inhibition_analysis_match_test$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

inhibition_model <- rma.mv(yi ~ comparison,
                           V = variance_matrix_inhibition,
                           random = ~ 1 | StudyID/ESID,
                           #subset = (inhibition_analysis_match_test$StudyID != 1), #for leave-1-out
                           data=inhibition_analysis_match_test)
inhibition_model_robust <- robust(inhibition_model, cluster = inhibition_analysis_match_test$StudyID, clubSandwich = TRUE)
predict(inhibition_model_robust, transf=exp, digits=2)
