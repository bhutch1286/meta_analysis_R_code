#----------Sensitivity analysis with low precision trimmed dataset

IB_summary <- summary(IB_dataset$sinv)
PB_lp_trimmed <- IB_dataset[!IB_dataset$sinv >= IB_summary[5],] #upper quartile
sum(IB_dataset$sinv >= IB_summary[5])
IB_dataset$trimmed <- IB_dataset$sinv >= IB_summary[5]

#recompute varcov matrix
variance_matrix_pboor <- 
  impute_covariance_matrix(PB_lp_trimmed$vi, 
                           cluster = PB_lp_trimmed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

#Checking if significant moderators (online, inhibition, targets, distractors) still show significance
ml_model_pubs_online2 <- rma.mv(yi ~ sinv + online,
                                V = variance_matrix_pboor,
                                random = ~ 1 | StudyID/ESID,
                                data=PB_lp_trimmed,
                                subset = (PB_lp_trimmed$literature == "AS"))
robust(ml_model_pubs_online2, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#inhibition
inhib_v_no_inhib <- rma.mv(yi ~ sinv + inhibition,
                           V = variance_matrix_pboor,
                           random = ~ 1 | StudyID/ESID,
                           data=PB_lp_trimmed,
                           subset = (PB_lp_trimmed$literature == "AS"))
robust(inhib_v_no_inhib, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)


#targets
ml_model_target_pubs_out_remo <- rma.mv(yi ~ sinv + targets,
                                        V=variance_matrix_pboor,
                                        random = ~ 1 | StudyID/ESID,
                                        data=PB_lp_trimmed,
                                        subset = (PB_lp_trimmed$literature == "AS"))
robust(ml_model_target_pubs_out_remo, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#distractors
ml_model_distractors_pubs_out_remo <- rma.mv(yi ~ sinv + distractors,
                                            V=variance_matrix_pboor,
                                            random = ~ 1 | StudyID/ESID,
                                            data=PB_lp_trimmed,
                                            subset = (PB_lp_trimmed$literature == "AS"))
robust(ml_model_distractors_pubs_out_remo, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)


#sinv is no longer significant in all >>>> proceed with remainder of sensitivity analysis


#main model
SA_pb_test <- rma.mv(yi,
                     V = variance_matrix_pboor,
                     random = ~ 1 | StudyID/ESID,
                     data=PB_lp_trimmed)
SA_pb_test_robust <- robust(SA_pb_test, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)
predict(SA_pb_test_robust, transf=exp, digits=2)

#literature
SA_lit_ml_model <- rma.mv(yi ~ literature,
                          V = variance_matrix_pboor,
                          random = ~ 1 | StudyID/ESID, 
                          data=PB_lp_trimmed)
SA_lit_ml_model_robust <- robust(SA_lit_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)
predict(SA_lit_ml_model_robust, transf=exp, digits=2)


#age
SA_age_ml_model <- rma.mv(yi ~ age_imputed,
                          V = variance_matrix_pboor,
                          random = ~ 1 | StudyID/ESID ,
                          data=PB_lp_trimmed)
robust(SA_age_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#gender
SA_gender_ml_model <- rma.mv(yi ~ male_prop_imputed,
                             V = variance_matrix_pboor,
                             random = ~ 1 | StudyID/ESID , 
                             data=PB_lp_trimmed)
robust(SA_gender_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#paradigm
SA_paradigm_ml_model <- rma.mv(yi ~ paradigm,
                               V = variance_matrix_pboor,
                               random = ~ 1 | StudyID/ESID ,
                               data=PB_lp_trimmed)
robust(SA_paradigm_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#literature paradigm interaction (for additive model, replace * with +)
SA_paradigm_ml_model <- rma.mv(yi ~ literature + paradigm,
                               V = variance_matrix_pboor,
                               random = ~ 1 | StudyID/ESID , 
                               data=PB_lp_trimmed)
robust(SA_paradigm_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#setting
SA_setting_ml_model <- rma.mv(yi ~ setting,
                              V = variance_matrix_pboor,
                              random = ~ 1 | StudyID/ESID , 
                              data=PB_lp_trimmed)
robust(SA_setting_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#task
SA_task_ml_model <- rma.mv(yi ~ task_imputed,
                           V = variance_matrix_pboor,
                           random = ~ 1 | StudyID/ESID , 
                           data=PB_lp_trimmed)
robust(SA_task_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#CS duration
PB_lp_trimmed_CS_d_outlier_removed <- PB_lp_trimmed[!PB_lp_trimmed$CS_d >= 15,]

variance_matrix_pbcsdor <- 
  impute_covariance_matrix(PB_lp_trimmed_CS_d_outlier_removed$vi, 
                           cluster = PB_lp_trimmed_CS_d_outlier_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

ml_model_pb_CSd_outliers_removed <- rma.mv(yi ~ CS_d,
                                           V = variance_matrix_pbcsdor,
                                           random = ~ 1 | Study/ESID,
                                           data=PB_lp_trimmed_CS_d_outlier_removed)
robust(ml_model_pb_CSd_outliers_removed, cluster=PB_lp_trimmed_CS_d_outlier_removed$StudyID, clubSandwich = TRUE)


#trial duration
PB_lp_trimmed_triald_outliers_removed <- PB_lp_trimmed[!PB_lp_trimmed$trial_d > 40,]

variance_matrix_pbtdor <- 
  impute_covariance_matrix(PB_lp_trimmed_triald_outliers_removed$vi, 
                           cluster = PB_lp_trimmed_triald_outliers_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

ml_model_pb_TD_outliers_removed <- rma.mv(yi ~ trial_d,
                                          V = variance_matrix_pbtdor,
                                          random = ~ 1 | StudyID/ESID,
                                          data=PB_lp_trimmed_triald_outliers_removed)
robust(ml_model_pb_TD_outliers_removed, cluster = PB_lp_trimmed_triald_outliers_removed$StudyID, clubSandwich = TRUE)


#CS proportion
PB_lp_trimmed_CSp_outlier_removed <- PB_lp_trimmed[!PB_lp_trimmed$CS_prop > 0.80,]

variance_matrix_pbcspor <- 
  impute_covariance_matrix(PB_lp_trimmed_CSp_outlier_removed$vi, 
                           cluster = PB_lp_trimmed_CSp_outlier_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

ml_model_CSp_outliers_removed <- rma.mv(yi ~ CS_prop,
                                        V = variance_matrix_pbcspor,
                                        random = ~ 1 | Study/ESID,
                                        data=PB_lp_trimmed_CSp_outlier_removed)
robust(ml_model_CSp_outliers_removed, cluster=PB_lp_trimmed_CSp_outlier_removed$StudyID, clubSandwich = TRUE)


#pre critical trials
PB_lp_trimmed$pre_ct_outliers <- PB_lp_trimmed$pre_ct > 15
PB_lp_prect_outliers <- which(PB_lp_trimmed$pre_ct_outliers)
PB_lp_prect_outliers_removed <- PB_lp_trimmed[-c(PB_lp_prect_outliers),]

variance_matrix_pbprector <- impute_covariance_matrix(PB_lp_prect_outliers_removed$vi, 
                                                    cluster = PB_lp_prect_outliers_removed$StudyID, 
                                                    r = rho,
                                                    smooth_vi = TRUE, check_PD = TRUE)

SA_preCT_outliers_rem_ml_model <- rma.mv(yi ~ pre_ct,
                                         V = variance_matrix_pbprector,
                                         random = ~ 1 | StudyID/ESID , 
                                         data = PB_lp_prect_outliers_removed)
robust(SA_preCT_outliers_rem_ml_model, cluster = PB_lp_prect_outliers_removed$StudyID, clubSandwich = TRUE)

#targets
SA_targets_ml_model <- rma.mv(yi ~ targets,
                              V = variance_matrix_pboor,
                              random = ~ 1 | StudyID/ESID , 
                              data=PB_lp_trimmed,
                              subset = (PB_lp_trimmed$literature == "AS"))
robust(SA_targets_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#distractors
SA_distractors_ml_model <- rma.mv(yi ~ distractors,
                                 V = variance_matrix_pboor,
                                 random = ~ 1 | StudyID/ESID , 
                                 data=PB_lp_trimmed,
                                 subset = (PB_lp_trimmed$literature == "AS"))
robust(SA_distractors_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#prior knowledge
SA_PK_ml_model <- rma.mv(yi ~ PK,
                         V = variance_matrix_pboor,
                         random = ~ 1 | StudyID/ESID , 
                         data=PB_lp_trimmed)
robust(SA_PK_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#FA trial
SA_FA_ml_model <- rma.mv(yi ~ FA,
                         V = variance_matrix_pboor,
                         random = ~ 1 | StudyID/ESID , 
                         data=PB_lp_trimmed)
robust(SA_FA_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#FA trial load only
SA_FA_load_ml_model <- rma.mv(yi ~ FA,
                              V = variance_matrix_pboor,
                              random = ~ 1 | StudyID/ESID , 
                              data=PB_lp_trimmed,
                              subset = (PB_lp_trimmed$literature == "load"))
robust(SA_FA_load_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#task performance
SA_TP_ml_model <- rma.mv(yi ~ TP,
                         V = variance_matrix_pboor,
                         random = ~ 1 | StudyID/ESID , 
                         data=PB_lp_trimmed)
robust(SA_TP_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#online
SA_online_ml_model <- rma.mv(yi ~ online, 
                             V = variance_matrix_pboor,
                             random = ~ 1 | StudyID/ESID , 
                             data=PB_lp_trimmed)
robust(SA_online_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#load validity
SA_load_val_ml_model <- rma.mv(yi ~ load_val,
                               V = variance_matrix_pboor,
                               random = ~ 1 | StudyID/ESID , 
                               data=PB_lp_trimmed)
robust(SA_load_val_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#quality score
SA_qscore_ml_model <- rma.mv(yi ~ q_score,
                             V = variance_matrix_pboor,
                             random = ~ 1 | StudyID/ESID , 
                             data=PB_lp_trimmed)
robust(SA_qscore_ml_model, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)



#load difference score
PB_lp_trimmed$load_diff_outliers <- PB_lp_trimmed$target_distractor_l_diff > 5
PB_lp_load_diff_outliers <- which(PB_lp_trimmed$load_diff_outliers)
load_diff_PB_lp_removed <- PB_lp_trimmed[-c(PB_lp_load_diff_outliers),]

#need to recompute var cov matrix
variance_matrix_ldpboor <- 
  impute_covariance_matrix(load_diff_PB_lp_removed$vi, 
                           cluster = load_diff_PB_lp_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

SA_load_difference_analysis_outlier_remv <- rma.mv(yi ~ target_distractor_l_diff,
                                                   V = variance_matrix_ldpboor,
                                                   random = ~ 1 | StudyID/ESID , 
                                                   subset = (load_diff_PB_lp_removed$literature == "load" &
                                                               load_diff_PB_lp_removed$o_l_manip == 1 &
                                                               load_diff_PB_lp_removed$task == "search"),
                                                   data=load_diff_PB_lp_removed)
robust(SA_load_difference_analysis_outlier_remv, cluster = load_diff_PB_lp_removed$StudyID, clubSandwich = TRUE)


#Attention set level
AS_ml_model_level <- rma.mv(yi ~ level,
                            V = variance_matrix_pboor,
                            random = ~ 1 | StudyID/ESID,
                            data=PB_lp_trimmed,
                            subset = (PB_lp_trimmed$literature == "AS"))
robust(AS_ml_model_level, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#Attention set type
AS_ml_model_type <- rma.mv(yi ~ type,
                           V = variance_matrix_pboor,
                           random = ~ 1 | StudyID/ESID,
                           data=PB_lp_trimmed,
                           subset = (PB_lp_trimmed$literature == "AS"))
robust(AS_ml_model_type, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)


#Load type
load_ml_model_type <- rma.mv(yi ~ load_type,
                             V = variance_matrix_pboor,
                             random = ~ 1 | StudyID/ESID,
                             data=PB_lp_trimmed,
                             subset = (PB_lp_trimmed$literature == "load"))
robust(load_ml_model_type, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#Inhibition
inhib_v_no_inhib <- rma.mv(yi ~ inhibition,
                           V = variance_matrix_pboor,
                           random = ~ 1 | StudyID/ESID , 
                           data=PB_lp_trimmed,
                           subset = (PB_lp_trimmed$literature == "AS"))
robust(inhib_v_no_inhib, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#published vs unpublished
ml_model_published <- rma.mv(yi ~ published,
                             V=variance_matrix_pboor,
                             random = ~ 1 | StudyID/ESID,
                             data=PB_lp_trimmed)
robust(ml_model_published, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)

#Peters regression
SA_pb_test <- rma.mv(yi ~ sinv,
                     V = variance_matrix_pboor,
                     random = ~ 1 | StudyID/ESID,
                     data=PB_lp_trimmed)
SA_pb_test_robust <- robust(SA_pb_test, cluster = PB_lp_trimmed$StudyID, clubSandwich = TRUE)



#----other datasets needed for these analyses

#inhibition analysis
inhibition_analysis_match_test <- inhibition_analysis[inhibition_analysis$matchv1 & !inhibition_analysis$low_prec_removed,] #unhash for low precision study trim analysis

rho <- 0.6
variance_matrix_inhibition <- 
  impute_covariance_matrix(inhibition_analysis_match_test$vi, 
                           cluster = inhibition_analysis_match_test$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

inhibition_model <- rma.mv(yi ~ comparison,
                           V = variance_matrix_inhibition,
                           random = ~ 1 | StudyID/ESID,
                           data=inhibition_analysis_match_test)
inhibition_model_robust <- robust(inhibition_model, cluster = inhibition_analysis_match_test$StudyID, clubSandwich = TRUE)



#interaction analysis

#model 1
AS_interaction <- rma.mv(yi ~ factor(HL),
                         V = varmat_AS_interaction, 
                         random =  ~ 1 | studyID/ESID,
                         subset = (interaction_AS$low_prec_removed != 1), 
                         data =  interaction_AS)
AS_interaction_robust <- robust(AS_interaction, cluster = interaction_AS$studyID, clubSandwich = TRUE)


#model 2
load_interaction <- rma.mv(yi ~ factor(match),
                           V = varmat_load_interaction,
                           random =  ~ 1 | studyID/ESID,
                           subset = (interaction_load$low_prec_removed != 1),
                           data = interaction_load)
load_interaction_robust <- robust(load_interaction, cluster = interaction_load$studyID, clubSandwich = TRUE)
