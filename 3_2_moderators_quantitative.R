#----------Quantitative moderators


#1.CS duration
ml_model_CSd <- rma.mv(yi ~ CS_d,
                       V = variance_matrix,
                       random = ~ 1 | Study/ESID,
                       data=IB_dataset)
robust(ml_model_CSd, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_CSd, mod="CS_d", xlab="CS duration", pred=TRUE, ci=TRUE)

#outliers
IB_v3_CS_d_outlier_removed <- IB_dataset[!IB_dataset$CS_d >= 15,]
sum(IB_dataset$CS_d >= 15)

#recompute varcov matrix
variance_matrix_csdor <- 
  impute_covariance_matrix(IB_v3_CS_d_outlier_removed$vi, 
                           cluster = IB_v3_CS_d_outlier_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

ml_model_CSd_outliers_removed <- rma.mv(yi ~ CS_d,
                                        V = variance_matrix_csdor,
                                        random = ~ 1 | Study/ESID,
                                        data=IB_v3_CS_d_outlier_removed)
robust(ml_model_CSd_outliers_removed, cluster=IB_v3_CS_d_outlier_removed$StudyID, clubSandwich = TRUE)
regplot(ml_model_CSd_outliers_removed, mod="CS_d", xlab="CS duration", pred=TRUE, ci=TRUE)


#2.trial duration
ml_model_TD <- rma.mv(yi ~ trial_d,
                      V = variance_matrix,
                      random = ~ 1 | StudyID/ESID,
                      data=IB_dataset)
robust(ml_model_TD, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_TD, mod="trial_d", xlab="Trial duration", pred=TRUE, ci=TRUE,label=TRUE)

#outliers
IB_v3_triald_outliers_removed <- IB_dataset[!IB_dataset$trial_d > 40,]
sum(IB_dataset$trial_d > 40)

#recompute varcov matrix
variance_matrix_tdor <- 
  impute_covariance_matrix(IB_v3_triald_outliers_removed$vi, 
                           cluster = IB_v3_triald_outliers_removed$StudyID, 
                           r = rho, smooth_vi = TRUE, check_PD = TRUE)

ml_model_TD_outliers_removed <- rma.mv(yi ~ trial_d,
                                       V = variance_matrix_tdor,
                                       random = ~ 1 | StudyID/ESID,
                                       data=IB_v3_triald_outliers_removed)
robust(ml_model_TD_outliers_removed, cluster=IB_v3_triald_outliers_removed$StudyID, clubSandwich = TRUE)
regplot(ml_model_TD_outliers_removed, mod="trial_d", xlab="Trial duration", pred=TRUE, ci=TRUE)


#3.CS proportion
ml_model_CSp <- rma.mv(yi ~ CS_prop,
                       V = variance_matrix,
                       random = ~ 1 | StudyID/ESID,
                       data=IB_dataset)
robust(ml_model_CSp, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_CSp, mod="CS_prop", xlab="CS proportion", pred=TRUE, ci=TRUE)

#outliers
IB_dataset_CSp_outlier_removed <- IB_dataset[!IB_dataset$CS_prop > 0.80,]
sum(IB_dataset$CS_prop > 0.80)

#recompute varcov matrix
variance_matrix_cspor <- 
  impute_covariance_matrix(IB_dataset_CSp_outlier_removed$vi, 
                           cluster = IB_dataset_CSp_outlier_removed$StudyID, 
                           r = rho, smooth_vi = TRUE, check_PD = TRUE)

ml_model_CSp_outliers_removed <- rma.mv(yi ~ CS_prop,
                                        V = variance_matrix_cspor,
                                        random = ~ 1 | Study/ESID,
                                        data=IB_dataset_CSp_outlier_removed)
robust(ml_model_CSp_outliers_removed, cluster=IB_dataset_CSp_outlier_removed$StudyID, clubSandwich = TRUE)
regplot(ml_model_CSp_outliers_removed, mod="CS_prop", xlab="CS proportion", ylab = "Log Odds", pred=TRUE, ci=TRUE)

#interaction with literature (interaction/additive model, replace * with +)
ml_model_CSp_outliers_removed_int <- rma.mv(yi ~ literature * CS_prop,
                                        V = variance_matrix_cspor,
                                        random = ~ 1 | Study/ESID,
                                        data=IB_dataset_CSp_outlier_removed)
robust(ml_model_CSp_outliers_removed_int, cluster=IB_dataset_CSp_outlier_removed$StudyID, clubSandwich = TRUE)

#AS literature
ml_model_CSp_outliers_removed_AS <- rma.mv(yi ~ CS_prop,
                                        V = variance_matrix_cspor,
                                        random = ~ 1 | Study/ESID,
                                        data=IB_dataset_CSp_outlier_removed,
                                        subset = (IB_dataset_CSp_outlier_removed$literature == "AS"))
ml_model_CSp_outliers_removed_AS_robust <- robust(ml_model_CSp_outliers_removed_AS, cluster=IB_dataset_CSp_outlier_removed$StudyID, clubSandwich = TRUE)
regplot(ml_model_CSp_outliers_removed_AS_robust, mod="CS_prop", xlab="CS proportion", ylab = "Log Odds", pred=TRUE, ci=TRUE)


#4. pre critical trials
ml_model_preCT <- rma.mv(yi ~ pre_ct,
                         V = variance_matrix,
                         random = ~ 1 | StudyID/ESID,
                         data=IB_dataset)
robust(ml_model_preCT, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(ml_model_preCT, mod="pre_ct", xlab="number of pre critical trials", pred=TRUE, ci=TRUE,label=TRUE)

#outliers
IB_dataset$pre_ct_outliers <- IB_dataset$pre_ct > 15
prect_outliers <- which(IB_dataset$pre_ct_outliers)
sum(IB_dataset$pre_ct > 15, na.rm = TRUE)

IB_v3_prect_outliers_removed <- IB_dataset[-c(prect_outliers ),]

#recompute varcov matrix
variance_matrix_prector <- impute_covariance_matrix(IB_v3_prect_outliers_removed$vi, 
                         cluster = IB_v3_prect_outliers_removed$StudyID, 
                         r = rho, smooth_vi = TRUE, check_PD = TRUE)

ml_model_preCT_outliers_removed <- rma.mv(yi ~  pre_ct,
                                          V = variance_matrix_prector,
                                          random = ~ 1 | StudyID/ESID,
                                          #subset = (IB_v3_prect_outliers_removed$literature == "AS"), #unhash for effect in attention set literature only
                                          data=IB_v3_prect_outliers_removed)
ml_model_preCT_outliers_removed_robust <- robust(ml_model_preCT_outliers_removed, cluster=IB_v3_prect_outliers_removed$StudyID, clubSandwich = TRUE)
regplot(ml_model_preCT_outliers_removed_robust, mod="pre_ct", xlab="number of pre critical trials", pred=TRUE, ci=TRUE)


#4.targets
AS_ml_model_target <- rma.mv(yi ~ targets,
                             V = variance_matrix,
                             random = ~ 1 | StudyID/ESID,
                             data=IB_dataset,
                             subset = (IB_dataset$literature == "AS"))
AS_ml_model_target_robust <- robust(AS_ml_model_target, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(AS_ml_model_target_robust, mod="target", xlab="Number of targets", ylab = "Log Odds", pred=TRUE, ci=TRUE)


#5.distractors
AS_ml_model_distractor <- rma.mv(yi ~ distractors,
                                 V = variance_matrix,
                                 random = ~ 1 | StudyID/ESID,
                                 data=IB_dataset,
                                 subset = (IB_dataset$literature == "AS"))
AS_ml_model_distractor_robust <- robust(AS_ml_model_distractor, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(AS_ml_model_distractor_robust, mod="distractor", xlab="Number of distractors", ylab = "Log Odds", pred=TRUE, ci=TRUE)
