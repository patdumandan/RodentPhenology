#Data Manipulation####
portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pb_plot=portal_reprod%>%filter(species=="PB")
PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")


pp_plot=portal_reprod%>%filter(species=="PP")
PP_male_con=pp_plot%>%filter(treatment=="control", sex=="male")
PP_male_ex=pp_plot%>%filter(treatment=="exclosure", sex=="male")
PP_female_con=pp_plot%>%filter(treatment=="control", sex=="female")
PP_female_ex=pp_plot%>%filter(treatment=="exclosure", sex=="female")

dipos_all=portal_reprod%>%filter(species%in%c("DO", "DM", "DS"))

dipo_f=dipos_all%>%filter(sex=="female")
dipo_m=dipos_all%>%filter(sex=="male")

#Treatment effect####
#data###
portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pbf_plot=portal_reprod%>%filter(sex=="female", species=="PB")
pbm_plot=portal_reprod%>%filter(sex=="male", species=="PB")
ppf_plot=portal_reprod%>%filter(sex=="female", species=="PP")
ppm_plot=portal_reprod%>%filter(sex=="male", species=="PP")

#build GAMs with treatment as ordered factors to assess treatment effect

pbf_plot=pbf_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
pbm_plot=pbm_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
ppf_plot=ppf_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
ppm_plot=ppm_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))

#PB models##

pbf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbf_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n",main="female", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

pbm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbm_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", main="male",ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#PP models##

ppf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppf_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

ppm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppm_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

summary(pbf3)
summary(pbm3)
summary(ppf3)
summary(ppm3)

#Abiotic and Biotic Drivers####

#PB models###
pbm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_con, method = 'REML', weights = abundance, family = binomial)
pbm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

pbf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_con, method = 'REML', weights = abundance, family = binomial)
pbf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

#PP models##
ppm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_con, method = 'REML', weights = abundance, family = binomial)
ppm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

ppf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_con, method = 'REML', weights = abundance, family = binomial)
ppf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

##k-rat models##
dipof_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipo_f, method = 'REML', weights = abundance, family = binomial)
dipom_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipo_m, method = 'REML', weights = abundance, family = binomial)


#results viz###
plot.gam(pbf3, shade=T, pages=1)
