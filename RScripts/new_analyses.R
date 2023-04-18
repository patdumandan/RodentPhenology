#data###
portal_reprod=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

pbf_plot=portal_reprod%>%filter(sex=="female", species=="PB")
pbm_plot=portal_reprod%>%filter(sex=="male", species=="PB")
ppf_plot=portal_reprod%>%filter(sex=="female", species=="PP")
ppm_plot=portal_reprod%>%filter(sex=="male", species=="PP")

pbf_plot=pbf_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
pbm_plot=pbm_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
ppf_plot=ppf_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))
ppm_plot=ppm_plot%>%mutate(otrt=ordered(treatment, levels=c("control", "exclosure")))

#PB models####

pbf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 ndvis+ndvis_lag+temps_mean+
                 temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ 
                 pbs+pps+dipos+otrt, data=pbf_plot, method = 'REML', weights = abundance, family = binomial)

pbm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 ndvis+ndvis_lag+temps_mean+
                 temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ 
                 pbs+pps+dipos+otrt, data=pbm_plot, method = 'REML', weights = abundance, family = binomial)

#PP models####

ppf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 ndvis+ndvis_lag+temps_mean+
                 temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ 
                 pbs+pps+dipos+otrt, data=ppf_plot, method = 'REML', weights = abundance, family = binomial)

ppm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 ndvis+ndvis_lag+temps_mean+
                 temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ 
                 pbs+pps+dipos+otrt, data=ppm_plot, method = 'REML', weights = abundance, family = binomial)

summary(pbf3)
summary(pbm3)
summary(ppf3)
summary(ppm3)

#results viz####
plot.gam(pbf3, shade=T, pages=1)
