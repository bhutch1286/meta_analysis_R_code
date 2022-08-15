sorted_IB <- IB_dataset[order(IB_dataset$YoP),]

variance_matrix <- 
  impute_covariance_matrix(sorted_IB$vi, 
                           cluster = sorted_IB$StudyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

#--------------------------------Load
forest_cog <- rma.mv(yi, V = variance_matrix,
                       random = ~ 1 | StudyID/ESID,
                       data=sorted_IB,
                       subset = (sorted_IB$load_type == "cognitive"))
forest_cog_robust <- robust(forest_cog, cluster = sorted_IB$StudyID, clubSandwich = TRUE)

forest(forest_cog_robust, addpre=TRUE, xlim=c(-8,4.5),
       slab = sorted_IB$Study, ilab = sorted_IB$Subgroup_2, ilab.xpos=c(-5.2),
       mlab = "Multilevel model with RVE", header="Author(s) and Year", xlab="Log Odds Ratio" )
text(c(-5.2), 42, c("Subgroup"), cex=.75, font=2)

forest_percep <- rma.mv(yi, V = variance_matrix,
                                random = ~ 1 | StudyID/ESID, data=sorted_IB,
                                subset = (sorted_IB$load_type == "perceptual"))
forest_percep_robust <- robust(forest_percep, cluster=sorted_IB$StudyID, clubSandwich = TRUE)

forest(forest_percep_robust, addpre=TRUE, xlim=c(-8,7),
       slab = sorted_IB$Study, ilab = sorted_IB$Subgroup_2, ilab.xpos=c(-5.2),
       mlab = "Multilevel model with RVE", header="Author(s) and Year", xlab="Log Odds Ratio" )
text(c(-5.2), 47, c("Subgroup"), cex=.6, font=2)


#--------------------------------Attention set

forest_feature <- rma.mv(yi, V = variance_matrix,
                     random = ~ 1 | StudyID/ESID,
                     data=sorted_IB,
                     subset = (sorted_IB$level == "feature"))
forest_feature_robust <- robust(forest_feature, cluster = sorted_IB$StudyID, clubSandwich = TRUE)

forest(forest_feature_robust, addpre=TRUE, xlim=c(-8,4.5),
       slab = sorted_IB$Study, ilab = sorted_IB$Subgroup_2, ilab.xpos=c(-5.2),
       mlab = "Multilevel model with RVE", header="Author(s) and Year", xlab="Log Odds Ratio")
text(c(-5.2), 97, c("Subgroup"), cex=.25, font=2)


forest_semantic <- rma.mv(yi, V = variance_matrix,
                        random = ~ 1 | StudyID/ESID, data=sorted_IB,
                        subset = (sorted_IB$level == "semantic"))
forest_semantic_robust <- robust(forest_semantic, cluster=sorted_IB$StudyID, clubSandwich = TRUE)

forest(forest_semantic_robust, addpre=TRUE, xlim=c(-8,4.5),
       slab = sorted_IB$Study, ilab = sorted_IB$Subgroup_2, ilab.xpos=c(-5.2),
       mlab = "Multilevel model with RVE", header="Author(s) and Year", xlab="Log Odds Ratio" )
text(c(-5.2), 67, c("Subgroup"), cex=.4, font=2)


forest_inherent <- rma.mv(yi, V = variance_matrix,
                        random = ~ 1 | StudyID/ESID, data=sorted_IB,
                        subset = (sorted_IB$level == "inherent"))
forest_inherent_robust <- robust(forest_inherent, cluster=sorted_IB$StudyID, clubSandwich = TRUE)

forest(forest_inherent_robust, addpre=TRUE, xlim=c(-8,4.5),
       slab = sorted_IB$Study, ilab = sorted_IB$Subgroup_2, ilab.xpos=c(-5.2),
       mlab = "Multilevel model with RVE", header="Author(s) and Year", xlab="Log Odds Ratio" )
text(c(-5.2), 70, c("Subgroup"), cex=.4, font=2)
