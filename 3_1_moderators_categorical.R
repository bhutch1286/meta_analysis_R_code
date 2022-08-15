#----------Categorical moderators


#1. paradigm
table(IB_dataset$literature, IB_dataset$paradigm)
chisq.test(IB_dataset$literature, IB_dataset$paradigm, correct = FALSE)

ml_model_paradigm <- rma.mv(yi ~ paradigm, 
                            V = variance_matrix,
                            random = ~ 1 | StudyID/ESID, 
                            data=IB_dataset)
robust(ml_model_paradigm, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#interaction with literature (for additive model, replace * with +)
ml_model_paradigm_lit <- rma.mv(yi ~ literature * paradigm,
                                V = variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset)
robust(ml_model_paradigm_lit, cluster=IB_dataset$StudyID, clubSandwich = TRUE)


#2.setting
ml_model_setting <- rma.mv(yi ~ setting,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset)
robust(ml_model_setting, cluster=IB_dataset$StudyID, clubSandwich = TRUE)


#3. task
ml_model_task_imputed <- rma.mv(yi ~ task_imputed, 
                                V = variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset)
ml_model_task_imputed_robust <- robust(ml_model_task_imputed, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_task_imputed_robust, transf=exp, digits=2)


#interaction with literature (for additive model, replace * with +)
ml_model_task_imputed <- rma.mv(yi ~ literature * task_imputed,
                                V = variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset)
ml_model_task_imputed_robust <- robust(ml_model_task_imputed, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
ml_model_task_imputed_robust

#--pairwise comparisons
#recalculate same model without intercept term to use as input for subgroup analysis below
ml_model_task_imputed <- rma.mv(yi ~ 0 + task_imputed, 
                                V = variance_matrix,
                                random = ~ 1 | StudyID/ESID,
                                data=IB_dataset)
ml_model_task_imputed_robust <- robust(ml_model_task_imputed, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

# Subgroups
count_v_discrim <- anova(ml_model_task_imputed_robust, L=c(1, -1, 0, 0))
count_v_other <- anova(ml_model_task_imputed_robust, L=c(1, 0, -1, 0))
count_v_search <- anova(ml_model_task_imputed_robust, L=c(1, 0, 0, -1))
search_v_other <- anova(ml_model_task_imputed_robust, L=c(0, 0, 1, -1))
search_v_discrim <- anova(ml_model_task_imputed_robust, L=c(0, 1, 0, -1))
discrim_v_other <- anova(ml_model_task_imputed_robust, L=c(0, 1, -1, 0))

# p-values bonferonni adjusted
task_pvals <- c(count_v_discrim$pval, count_v_other$pval, count_v_search$pval, search_v_other$pval, search_v_discrim$pval, discrim_v_other$pval)
p.adjust(task_pvals, method = "bonferroni", n = length(task_pvals))



#----------Study quality moderators

#4. prior knowledge
ml_model_PK <- rma.mv(yi ~ PK,
                      V = variance_matrix,
                      random = ~ 1 | StudyID/ESID,
                      data=IB_dataset)
robust(ml_model_PK, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#5. Full attention trial
ml_model_FAt <- rma.mv(yi ~ FA,
                       V = variance_matrix,
                       random = ~ 1 | StudyID/ESID,
                       data=IB_dataset)
robust(ml_model_FAt, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#AS literature
ml_model_FAt_AS <- rma.mv(yi ~  FA,
                          V = variance_matrix,
                          random = ~ 1 | StudyID/ESID, data=IB_dataset,
                          subset = (IB_dataset$literature == "AS"))
robust(ml_model_FAt_AS, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#load literature
ml_model_FAt_load <- rma.mv(yi ~ FA,
                            V = variance_matrix,
                            random = ~ 1 | StudyID/ESID, data=IB_dataset,
                            subset = (IB_dataset$literature == "load"))
robust(ml_model_FAt_load, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_FAt_load, transf=exp, digits=2)

#6. Task performance
ml_model_TP <- rma.mv(yi ~ TP,
                      V = variance_matrix,
                      random = ~ 1 | StudyID/ESID,
                      data=IB_dataset)
robust(ml_model_TP, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#7. Online study (AS only)
ml_model_online <- rma.mv(yi ~ online,
                          V = variance_matrix,
                          random = ~ 1 | StudyID/ESID,
                          data=IB_dataset,
                          subset = (IB_dataset$literature == "AS"))
robust(ml_model_online, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#8. Load validity (load only)
ml_model_load_val <- rma.mv(yi ~ load_val,
                            V = variance_matrix,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_dataset)
robust(ml_model_load_val, cluster=IB_dataset$StudyID, clubSandwich = TRUE)


#9. Quality score
ml_model_q_score <- rma.mv(yi ~ q_score,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset)
robust(ml_model_q_score, cluster=IB_dataset$StudyID, clubSandwich = TRUE)

#exploratory analysis: q-score 100% vs q-score 0%
ml_model_q_score <- rma.mv(yi ~ q_score,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset,
                           subset = (IB_dataset$q_score == "100%_q" | IB_dataset$q_score == "0%_q"))
ml_model_q_score_robust <- robust(ml_model_q_score, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_q_score_robust, transf=exp, digits =2)





#----------Moderators unique to the Attention Set literature


#10. Attention set type
ml_model_AS_type <- rma.mv(yi ~  type,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset,
                           subset = (IB_dataset$literature == "AS"))
ml_model_AS_type_robust <- robust(ml_model_AS_type, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_AS_type_robust, transf=exp, digits =2)


#11. Attention set level
ml_model_AS_level <- rma.mv(yi ~ level,
                            V = variance_matrix,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_dataset,
                            subset = (IB_dataset$literature == "AS"))
ml_model_AS_level_robust <- robust(ml_model_AS_level, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_AS_level_robust, transf=exp, digits =2)

#--pairwise comparisons
#recalculate same model without intercept term to use as input for subgroup analysis below
ml_model_AS_level <- rma.mv(yi ~ 0 + level,
                            V = variance_matrix,
                            random = ~ 1 | StudyID/ESID,
                            data=IB_dataset,
                            subset = (IB_dataset$literature == "AS"))
ml_model_AS_level_robust <- robust(ml_model_AS_level, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_AS_level_robust, transf=exp, digits =2)

# Subgroups
feature_v_inherent <- anova(ml_model_AS_level_robust, L=c(1, -1, 0))
feature_v_semantic <- anova(ml_model_AS_level_robust, L=c(1, 0, -1))
inherent_v_semantic <- anova(ml_model_AS_level_robust, L=c(0, 1, -1))

# p-values bonferonni adjusted
AS_level_pvals <- c(feature_v_inherent$pval, feature_v_semantic$pval, inherent_v_semantic$pval)
p.adjust(AS_level_pvals, method = "bonferroni", n = length(AS_level_pvals))


#----------Moderators unique to the Load literature

#12. Load type
ml_model_load_type <- rma.mv(yi ~ load_type,
                             V = variance_matrix,
                             random = ~ 1 | StudyID/ESID,
                             data=IB_dataset,
                             subset = (IB_dataset$literature == "load"))
ml_model_load_type_robust <- robust(ml_model_load_type, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_load_type_robust, transf=exp, digits=2)

#13. Load sensitivity analysis
ml_model_load_sensit <- rma.mv(yi ~ SA_load, V = variance_matrix,
                               random = ~ 1 | StudyID/ESID,
                               data=IB_dataset,
                               subset = (IB_dataset$literature == "load"))
ml_model_load_sensit_robust <- robust(ml_model_load_sensit, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_model_load_sensit_robust, transf=exp, digits=2)

#14. Load difference score [Not a categorical moderator]
load_difference_analysis <- rma.mv(yi ~ target_distractor_l_diff,
                                   V = variance_matrix,
                                   random = ~ 1 | StudyID/ESID , 
                                   subset = (IB_dataset$literature == "load" &
                                               IB_dataset$o_l_manip == 1 &
                                               IB_dataset$task == "search"),
                                   data=IB_dataset)
robust(load_difference_analysis, cluster=IB_dataset$StudyID, clubSandwich = TRUE)
regplot(load_difference_analysis, pred=TRUE, ci=TRUE)


#Removing outliers
IB_dataset$load_diff_outliers <- IB_dataset$target_distractor_l_diff > 5
load_diff_outliers <- which(IB_dataset$load_diff_outliers)
load_diff_outlier_removed <- IB_dataset[-c(load_diff_outliers),]

#recompute var cov matrix
variance_matrix_ldor <- 
  impute_covariance_matrix(load_diff_outlier_removed$vi, 
                           cluster = load_diff_outlier_removed$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

#Load difference score outliers removed
load_difference_analysis_outlier_removed <- rma.mv(yi ~ target_distractor_l_diff,
                                                   V = variance_matrix_ldor,
                                                   random = ~ 1 | StudyID/ESID , 
                                                   subset = (load_diff_outlier_removed$literature == "load" &
                                                             load_diff_outlier_removed$o_l_manip == 1 &
                                                             load_diff_outlier_removed$task == "search"),
                                                data=load_diff_outlier_removed)
robust(load_difference_analysis_outlier_removed, cluster=load_diff_outlier_removed$StudyID, clubSandwich = TRUE)
regplot(load_difference_analysis_outlier_removed, pred=TRUE, ci=TRUE)
