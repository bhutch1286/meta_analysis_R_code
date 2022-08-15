#----------Interaction analysis


#1. Attention set effect size, load subgroup

interaction_AS <- read_excel("C:/Documents/interaction_AS_by_load_subgroup.xlsx")

#-----reverse the subgroups of cognitive load studies, only unhash if running cog load reverse coded
#interaction_AS$HL[interaction_AS$studyID == 8] <- c(0,0,0,1,1,1,1,1,1)

rho <- 0.6
varmat_AS_interaction <- 
  impute_covariance_matrix(interaction_AS$vi, 
                           cluster = interaction_AS$studyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

AS_interaction <- rma.mv(yi ~ factor(HL),
                         V = varmat_AS_interaction, 
                         random =  ~ 1 | studyID/ESID,
                         #subset = (interaction_AS$low_prec_removed != 1), #unhash if running low precision trimmed analysis
                         data =  interaction_AS)
AS_interaction_robust <- robust(AS_interaction, cluster = interaction_AS$studyID, clubSandwich = TRUE)
predict(AS_interaction_robust, transf=exp, digits=2)

#no intercept model for significance of AS (p value) in low and high load subgroups separately
AS_interaction_noi <- rma.mv(yi ~ 0 + factor(HL),
                             V = varmat_AS_interaction,
                             random =  ~ 1 | studyID/ESID,
                             data = interaction_AS)
AS_interaction_noi_robust <- robust(AS_interaction_noi, cluster = interaction_AS$studyID, clubSandwich = TRUE)
AS_interaction_noi_robust


#2. Load effect size, attention set subgroup

interaction_load <- read_excel("C:/Documents/interaction_load_by_AS_subgroup.xlsx")

#-----reverse effect size direction for cognitive load studies, only unhash if running cog load reverse coded
#CL_studies <- interaction_load$studyID == 6
#CL_studies[is.na(CL_studies)] = FALSE
#interaction_load$yi[CL_studies] = -interaction_load$yi[CL_studies]

rho <- 0.6
varmat_load_interaction <- 
  impute_covariance_matrix(interaction_load$vi, 
                           cluster = interaction_load$studyID, 
                           r = rho,
                           smooth_vi = TRUE, check_PD = TRUE)

load_interaction <- rma.mv(yi ~ factor(match),
                           V = varmat_load_interaction,
                           random =  ~ 1 | studyID/ESID,
                           #subset = (interaction_load$low_prec_removed != 1), #unhash if running low precision trimmed analysis
                           data = interaction_load)
load_interaction_robust <- robust(load_interaction, cluster = interaction_load$studyID, clubSandwich = TRUE)
predict(load_interaction_robust, transf=exp, digits=2)


#no intercept model for significance (p value) of load in match and mismatch subgroups separately
load_interaction_noi <- rma.mv(yi ~ 0 + factor(match),
                               V = varmat_load_interaction,
                               random =  ~ 1 | studyID/ESID,
                               data = interaction_load)
load_interaction_noi_robust <- robust(load_interaction_noi, cluster = interaction_load$studyID, clubSandwich = TRUE)
load_interaction_noi_robust
