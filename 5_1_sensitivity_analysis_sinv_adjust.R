#----------Sensitivity analysis: Model re-estimates adjusting for inverse sample size

#main model
ml_peters_model <- rma.mv(yi ~ sinv,
                          V=variance_matrix,
                          random = ~ 1 | StudyID/ESID,
                          data=IB_dataset)
ml_peters_model_robust <- robust(ml_peters_model, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

exp(ml_peters_model_robust$beta[1]) #exp OR
exp(ml_peters_model_robust$ci.lb[1]) #exp CI lower bound
exp(ml_peters_model_robust$ci.ub[1]) #exp CI upper bound

#literature (attention set and load separately) estimates are taken from the literature moderator model
ml_peters_model_lit <- rma.mv(yi ~ 0 + sinv + literature,
                              V=variance_matrix,
                              random = ~ 1 | StudyID/ESID,
                              data=IB_dataset)
ml_peters_model_lit_robust <- robust(ml_peters_model_lit, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

exp(ml_peters_model_lit_robust$beta[2]) #exp OR
exp(ml_peters_model_lit_robust$ci.lb[2]) #exp CI lower bound
exp(ml_peters_model_lit_robust$ci.ub[2]) #exp CI upper bound

exp(ml_peters_model_lit_robust$beta[3]) #exp OR
exp(ml_peters_model_lit_robust$ci.lb[3]) #exp CI lower bound
exp(ml_peters_model_lit_robust$ci.ub[3]) #exp CI upper bound




#----------Assessing for simultaneous moderation and publication bias

#Literature
ml_model_pubs_lit <- rma.mv(yi ~ sinv + literature,
                            V=variance_matrix,
                            random = ~ 1 | StudyID/ESID, 
                            data=IB_dataset)
ml_model_pubs_lit_robust <- robust(ml_model_pubs_lit, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#age
ml_model_pubs_age <- rma.mv(yi ~ sinv + age_imputed,
                            V=variance_matrix,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_dataset)
robust(ml_model_pubs_age, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#gender
ml_model_pubs_gender <- rma.mv(yi ~ sinv + male_prop_imputed,
                               V=variance_matrix,
                               random = ~ 1 | StudyID/ESID,
                               data=IB_dataset)
robust(ml_model_pubs_gender, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#paradigm
ml_model_pubs_paradigm <- rma.mv(yi ~ sinv + paradigm,
                                 V=variance_matrix,
                                 random = ~ 1 | StudyID/ESID,
                                 data=IB_dataset)
robust(ml_model_pubs_paradigm, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#setting
ml_model_pubs_setting <- rma.mv(yi ~ sinv + setting,
                                V=variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset)
robust(ml_model_pubs_setting, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#task
ml_model_pubs_task_imputed <- rma.mv(yi ~ sinv + task_imputed,
                                     V=variance_matrix,
                                     random = ~ 1 | StudyID/ESID,
                                     data=IB_dataset)
robust(ml_model_pubs_task_imputed, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#prior knowledge
ml_model_pubs_PK <- rma.mv(yi ~ sinv + PK,
                           V=variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset)
robust(ml_model_pubs_PK, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#FA trial
ml_model_pubs_FAt <- rma.mv(yi ~ sinv + FA,
                            V=variance_matrix,
                            random = ~ 1 | StudyID/ESID, 
                            data=IB_dataset)
robust(ml_model_pubs_FAt, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#FA trial load only
ml_model_pubs_FAt <- rma.mv(yi ~ sinv + FA,
                            V=variance_matrix,
                            random = ~ 1 | StudyID/ESID, 
                            data=IB_dataset,
                            subset = (IB_dataset$literature == "load"))
robust(ml_model_pubs_FAt, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#task performance
ml_model_pubs_TP <- rma.mv(yi ~ sinv + TP,
                           V=variance_matrix,
                           random = ~ 1 | StudyID/ESID, 
                           data=IB_dataset)
robust(ml_model_pubs_TP, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#quality score
ml_model_pubs_q_score <- rma.mv(yi ~ sinv + q_score,
                                V=variance_matrix,
                                random = ~ 1 | StudyID/ESID, 
                                data=IB_dataset)
robust(ml_model_pubs_q_score, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#load validity
ml_model_pubs_load_val <- rma.mv(yi ~ sinv + load_val,
                                 V=variance_matrix,
                                 random = ~ 1 | StudyID/ESID, 
                                 data=IB_dataset)
robust(ml_model_pubs_load_val, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#online
ml_model_pubs_online <- rma.mv(yi ~ sinv + online,
                               V=variance_matrix,
                               random = ~ 1 | StudyID/ESID,
                               data=IB_dataset,
                               subset = (IB_dataset$literature == "AS"))
robust(ml_model_pubs_online, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#CS duration
ml_model_pubs_CSd <- rma.mv(yi ~ sinv + CS_d,
                            V=variance_matrix_csdor,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_v3_CS_d_outlier_removed)
robust(ml_model_pubs_CSd, cluster = IB_v3_CS_d_outlier_removed$StudyID, clubSandwich = TRUE)

#trial duration
ml_model_pubs_triald <- rma.mv(yi ~ sinv + trial_d,
                               V=variance_matrix_tdor,
                               random = ~ 1 | StudyID/ESID,
                               data=IB_v3_triald_outliers_removed)
robust(ml_model_pubs_triald, cluster = IB_v3_triald_outliers_removed$StudyID, clubSandwich = TRUE)

#CS proportion
ml_model_pubs_CSp <- rma.mv(yi ~ sinv + CS_prop,
                            V=variance_matrix_cspor,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_dataset_CSp_outlier_removed)
robust(ml_model_pubs_CSp, cluster = IB_dataset_CSp_outlier_removed$StudyID, clubSandwich = TRUE)

#pre critical trials
ml_model_pubs_preCT_outliers_removed <- rma.mv(yi ~ sinv + pre_ct,
                                               V=variance_matrix_prector,
                                               random = ~ 1 | StudyID/ESID,
                                               data=IB_v3_prect_outliers_removed)
robust(ml_model_pubs_preCT_outliers_removed, cluster = IB_v3_prect_outliers_removed$StudyID, clubSandwich = TRUE)

#targets
ml_model_target_pubs_out_remo <- rma.mv(yi ~ sinv + targets,
                                        V=variance_matrix,
                                        random = ~ 1 | StudyID/ESID,
                                        data=IB_dataset,
                                        subset = (IB_dataset$literature == "AS"))
robust(ml_model_target_pubs_out_remo, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#distractors
ml_model_distractor_pubs_out_remo <- rma.mv(yi ~ sinv + distractors,
                                            V=variance_matrix,
                                            random = ~ 1 | StudyID/ESID,
                                            data=IB_dataset,
                                            subset = (IB_dataset$literature == "AS"))
robust(ml_model_distractor_pubs_out_remo, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#Attention set level
ml_model_pubs_AS_level <- rma.mv(yi ~ sinv + level,
                                 V=variance_matrix,
                                 random = ~ 1 | StudyID/ESID,
                                 data=IB_dataset,
                                 subset = (IB_dataset$literature == "AS"))
robust(ml_model_pubs_AS_level, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#Attention set type
ml_model_pubs_AS_type <- rma.mv(yi ~ sinv + type,
                                V=variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset,
                                subset = (IB_dataset$literature == "AS"))
robust(ml_model_pubs_AS_type, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#inhibition
inhib_v_no_inhib <- rma.mv(yi ~ sinv + inhibition,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset)
inhib_robust <- robust(inhib_v_no_inhib, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#load type
load_ml_model_type <- rma.mv(yi ~  sinv + load_type,
                             V = variance_matrix,
                             random = ~ 1 | StudyID/ESID, data=IB_dataset,
                             subset = (IB_dataset$literature == "load"))
load_ml_model_type_robust <- robust(load_ml_model_type, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#load difference score
load_difference_analysis_pubs_outlier_remv <- rma.mv(yi ~ sinv + target_distractor_l_diff,
                                                     V=variance_matrix_ldor,
                                                     random = ~ 1 | StudyID/ESID , 
                                                     subset = (load_diff_outlier_removed$literature == "load" &
                                                                 load_diff_outlier_removed$o_l_manip == 1 &
                                                                 load_diff_outlier_removed$task == "search"),
                                                     data=load_diff_outlier_removed)
robust(load_difference_analysis_pubs_outlier_remv, cluster = load_diff_outlier_removed$StudyID, clubSandwich = TRUE)

#published vs unpublished
ml_model_published <- rma.mv(yi ~ sinv + published,
                             V=variance_matrix,
                             random = ~ 1 | StudyID/ESID,
                             data=IB_dataset)
robust(ml_model_published, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
