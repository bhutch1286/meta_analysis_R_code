#----------Sensitivity analysis: worst case meta-analysis

Attention_set_dataset <- IB_dataset[!IB_dataset$literature == "load",]
Load_dataset <- IB_dataset[!IB_dataset$literature == "AS",]

S.est.whole.null <- svalue(yi = IB_dataset$yi,
                    vi = IB_dataset$vi,
                    q = log(1), model = "robust", favor.positive = TRUE,
                    clustervar = IB_dataset$StudyID,
                    eta.grid.hi = 200, return.worst.meta = TRUE)

exp(S.est.whole.null$meta.worst$reg_table$b.r) #exp OR
exp(S.est.whole.null$meta.worst$reg_table$CI.L) #exp CI lower bound
exp(S.est.whole.null$meta.worst$reg_table$CI.U) #exp CI upper bound

S.est.AS.null <- svalue(yi = Attention_set_dataset$yi,
                    vi = Attention_set_dataset$vi,
                    q = log(1), model = "robust", favor.positive = TRUE,
                    clustervar = Attention_set_dataset$StudyID,
                    eta.grid.hi = 200, return.worst.meta = TRUE)

exp(S.est.AS.null$meta.worst$reg_table$b.r) #exp OR
exp(S.est.AS.null$meta.worst$reg_table$CI.L) #exp CI lower bound
exp(S.est.AS.null$meta.worst$reg_table$CI.U) #exp CI upper bound

S.est.load.null <- svalue(yi = Load_dataset$yi,
                    vi = Load_dataset$vi,
                    q = log(1), model = "robust", favor.positive = TRUE,
                    clustervar = Load_dataset$StudyID,
                    eta.grid.hi = 200, return.worst.meta = TRUE)

exp(S.est.load.null$meta.worst$reg_table$b.r) #exp OR
exp(S.est.load.null$meta.worst$reg_table$CI.L) #exp CI lower bound
exp(S.est.load.null$meta.worst$reg_table$CI.U) #exp CI upper bound
