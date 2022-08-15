#----------Publication bias

#1.Funnel plots

#first need a no intercept model for funnel plot by literature ref line
ml_model_lit_noi <- rma.mv(yi ~ 0 + literature,
                           V = variance_matrix,
                           random = ~ 1 | StudyID/ESID,
                           data=IB_dataset)

funnel(x = IB_dataset$yi,
       ni=IB_dataset$sample,
       yaxis = "ninv",
       refline = ml_model$b,
       digits = 3L, xlab = "Log Odds")

funnel(x = IB_dataset$yi,
       ni=IB_dataset$sample,
       subset = (IB_dataset$literature == "AS"),
       yaxis = "ninv",
       refline = ml_model_lit_noi$b[1],
       digits = 3L, xlab = "Log Odds")

funnel(x = IB_dataset$yi,
       ni=IB_dataset$sample,
       subset = (IB_dataset$literature == "load"),
       yaxis = "ninv",
       refline = ml_model_lit_noi$b[2],
       digits = 3L, xlab = "Log Odds")


#2. published vs unpublished:
ml_model_published <- rma.mv(yi ~ published,
                             V=variance_matrix,
                             random = ~ 1 | StudyID/ESID,
                             data=IB_dataset)
robust(ml_model_published, cluster = IB_dataset$StudyID, clubSandwich = TRUE)


#3. Peter's regression
ml_peters_model <- rma.mv(yi ~ sinv,
                          V=variance_matrix,
                          random = ~ 1 | StudyID/ESID,
                          data=IB_dataset)
ml_peters_model_robust <- robust(ml_peters_model, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
predict(ml_peters_model_robust, transf=exp, digits = 2)

#Attention set literature only
ml_peters_model_AS <- rma.mv(yi ~ sinv,
                             V=variance_matrix,
                             random = ~ 1 | StudyID/ESID ,
                             data=IB_dataset,
                             subset = (IB_dataset$literature == "AS"))
robust(ml_peters_model_AS, cluster = IB_dataset$StudyID, clubSandwich = TRUE)

#Load literature only
ml_peters_model_load <- rma.mv(yi ~ sinv, V=variance_matrix,
                               random = ~ 1 | StudyID/ESID ,
                               data=IB_dataset,
                               subset = (IB_dataset$literature == "load"))
robust(ml_peters_model_load, cluster = IB_dataset$StudyID, clubSandwich = TRUE)
