#PB males####
PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

#control####

#abiotic only models####

#all variables with lag of 0 and lag of 1
pbcon_m1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                    temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                  data=PB_male_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
pbcon_m1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                     ppts_warm+ppts_cool,
                   data=PB_male_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
pbcon_m1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                     temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                   data=PB_male_con, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
pbcon_m1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=PB_male_con, weights=abundance, family=binomial(link="logit"))

summary(pbcon_m1_1) #13.4%
summary(pbcon_m1_2) #15.1%
summary(pbcon_m1_3) #14.5%
summary(pbcon_m1_4) #17.1%

#abiotic + competition models#####

pbcon_m2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2+pbs,
                     data=PB_male_con, weights=abundance, family=binomial(link="logit"))
summary(pbcon_m2) # 16.6%

pbcon_m3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2+pbs+pps+dipos,
                     data=PB_male_con, weights=abundance, family=binomial(link="logit"))

summary(pbcon_m3) # 20.2%

#exclosure####

#abiotic only models####

#all variables with lag of 0 and lag of 1
pbex_m1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                     data=PB_male_ex, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
pbex_m1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                       ppts_warm+ppts_cool,
                     data=PB_male_ex, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
pbex_m1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                     data=PB_male_ex, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
pbex_m1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=PB_male_ex, weights=abundance, family=binomial(link="logit"))

summary(pbex_m1_1) #35.2%
summary(pbex_m1_2) #32.7%
summary(pbex_m1_3) #35.1%
summary(pbex_m1_4) #38.8%

#abiotic + competition models#####

pbex_m2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                     temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2+pbs,
                   data=PB_male_ex, weights=abundance, family=binomial(link="logit"))
summary(pbex_m2) # 40.1%

pbex_m3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                     temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2+pbs+pps+dipos,
                   data=PB_male_ex, weights=abundance, family=binomial(link="logit"))

summary(pbex_m3) # 41.7%

#PB females####

#abiotic only models####

#all variables with lag of 0 and lag of 1
pbcon_f1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
pbcon_f1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                       ppts_warm+ppts_cool,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
pbcon_f1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
pbcon_f1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))

summary(pbcon_f1_1) #59%
summary(pbcon_f1_2) #53%
summary(pbcon_f1_3) #53.5%
summary(pbcon_f1_4) #50.2%

#abiotic + competition models####
pbcon_f2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+pbs,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))
pbcon_f3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+pbs+pps+dipos,
                     data=PB_female_con, weights=abundance, family=binomial(link="logit"))

summary(pbcon_f2) #58.9%
summary(pbcon_f3)#62.8%

#compare with other species####

#PP males####
PP_male_con=pp_plot%>%filter(treatment=="control", sex=="male")
PP_male_ex=pp_plot%>%filter(treatment=="exclosure", sex=="male")
PP_female_con=pp_plot%>%filter(treatment=="control", sex=="female")
PP_female_ex=pp_plot%>%filter(treatment=="exclosure", sex=="female")

#control####

#abiotic only models####

#all variables with lag of 0 and lag of 1
ppcon_m1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                     data=PP_male_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
ppcon_m1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                       ppts_warm+ppts_cool,
                     data=PP_male_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
ppcon_m1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                     data=PP_male_con, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
ppcon_m1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=PP_male_con, weights=abundance, family=binomial(link="logit"))

summary(ppcon_m1_1) #63.7%
summary(ppcon_m1_2) #61.7%
summary(ppcon_m1_3) #60.8%
summary(ppcon_m1_4) #58.8%

#abiotic + competition models#####

ppcon_m2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+pps,
                     data=PP_male_con, weights=abundance, family=binomial(link="logit"))
summary(ppcon_m2) # 64%

ppcon_m3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                     temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+pps+pbs+dipos,
                   data=PP_male_con, weights=abundance, family=binomial(link="logit"))
summary(ppcon_m3) # 70.4%

#PP females####
#exclosure####

#abiotic only models####

#all variables with lag of 0 and lag of 1
ppex_f1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                     data=PP_female_ex, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
ppex_f1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                       ppts_warm+ppts_cool,
                     data=PP_female_ex, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
ppex_f1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                     data=PP_female_ex, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
ppex_f1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=PP_female_ex, weights=abundance, family=binomial(link="logit"))

summary(ppex_f1_1) 
summary(ppex_f1_2) 
summary(ppex_f1_3) 
summary(ppex_f1_4) 


#DM females####

#control####

DM_male_con=dm_plot%>%filter(treatment=="control", sex=="male")
DM_male_ex=dm_plot%>%filter(treatment=="exclosure", sex=="male")
DM_female_con=dm_plot%>%filter(treatment=="control", sex=="female")
DM_female_ex=dm_plot%>%filter(treatment=="exclosure", sex=="female")

#all variables with lag of 0 and lag of 1
dmcon_f1_1=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+
                       temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool,
                     data=DM_female_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 0
dmcon_f1_2=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis+temps_mean+
                       ppts_warm+ppts_cool,
                     data=DM_female_con, weights=abundance, family=binomial(link="logit"))

#all variables with lag of 1
dmcon_f1_3=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm+ppts_lag_cool,
                     data=DM_female_con, weights=abundance, family=binomial(link="logit"))

#temp and NDVI with lag of 1 and precip with lag of 2
dmcon_f1_4=mgcv::gam(proportion~s(month, bs="cc")+s(year)+ndvis_lag+
                       temps_lag_mean+ppts_lag_warm_2+ppts_lag_cool_2,
                     data=DM_female_con, weights=abundance, family=binomial(link="logit"))

summary(dmcon_f1_1) #31.8%
summary(dmcon_f1_2) #29.5%
summary(dmcon_f1_3) #53.5%
summary(dmcon_f1_4) #50.2%
