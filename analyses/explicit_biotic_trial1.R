#PB MALES####
#control####
pb_plot=all3_bmass%>%filter(species=="PB")
pp_plot=all3_bmass%>%filter(species=="PP")
dm_plot=all3_bmass%>%filter(species=="DM")

PB_male_con=pb_plot%>%filter(treatment=="control", sex=="male")
PB_male_ex=pb_plot%>%filter(treatment=="exclosure", sex=="male")
PB_female_con=pb_plot%>%filter(treatment=="control", sex=="female")
PB_female_ex=pb_plot%>%filter(treatment=="exclosure", sex=="female")

#no explicit biotic effects
m1=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
                   , data=PB_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m1) #R2=40%

#add effect of DM biomass
m1_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)
                   , data=PB_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m1_bmass) #R2=42.7%

#add effect of DM biomass AND intraspecific competition (PB)
m1_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                   s(bmass_PB), data=PB_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m1_bmass2) #R2=48.5%

#best fit plot####
plot(m1_bmass2, shade=T, pages=1)
want=seq(1, nrow(PB_male_con), length.out = 200)
pdat2=with(PB_male_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                   precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                   bmass_PB[want]))
p3=predict(m1_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m1_bmass)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(m1_bmass2)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat2$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term)
abline(v=3.5, type="l", lty=2)

#exclosure####
#no explicit biotic effects
m2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m2) #R2=28.9%

#add effect of DM biomass
m2_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)
                   , data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m2_bmass) #R2=30.9%

#add effect of DM biomass AND intraspecific competition (PB competition)
m2_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB), data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m2_bmass2) #R2=32%

#best fit plot####
plot(m2_bmass2, shade=T, pages=1)
want=seq(1, nrow(PB_male_ex), length.out = 200)
pdat2=with(PB_male_ex, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                   precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                   bmass_PB[want]))
p3=predict(m2_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m2_bmass2)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m2.d <- Deriv(m2_bmass2)
Term="month"
m2.dci <- confint(m2.d, term = Term)
m2.dsig <- signifD(pdat2$p3, d = m2.d[[Term]]$deriv,
                   +m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
plot.Deriv(m2.d, sizer=T, term=Term)
abline(v=3.4, type="l", lty=2)


#PB FEMALES####

#control####
#no explicit biotic effects
m3=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PB_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m3) #R2=39.9%

#add effect of DM biomass
m3_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)
                   , data=PB_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m3_bmass) #R2=44.5%

#add effect of DM biomass AND intraspecific competition (PB competition)
m3_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB), data=PB_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m3_bmass2) #R2=44.5%

#best fit plot####
plot(m3_bmass, shade=T, pages=1)
want=seq(1, nrow(PB_female_con), length.out = 200)
pdat2=with(PB_female_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                  precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                  bmass_PB[want]))
p3=predict(m3_bmass, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m3_bmass)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m3.d <- Deriv(m3_bmass)
Term="month"
m3.dci <- confint(m3.d, term = Term)
m3.dsig <- signifD(pdat2$p3, d = m3.d[[Term]]$deriv,
                   +m3.dci[[Term]]$upper, m3.dci[[Term]]$lower)
plot.Deriv(m3.d, sizer=T, term=Term)
abline(v=4.6, type="l", lty=2)

#exclosure####

#no explicit biotic effects
m4=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m4) #R2=53.4%

#add effect of DM biomass
m4_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)
                   , data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m4_bmass) #R2=58.3%

#add effect of DM biomass AND intraspecific competition (PB competition)
m4_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB), data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m4_bmass2) #R2=58.5%

#best fit plot####
plot(m4_bmass, shade=T, pages=1)
want=seq(1, nrow(PB_female_ex), length.out = 200)
pdat2=with(PB_female_ex, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                  precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                  bmass_PB[want]))
p3=predict(m4_bmass, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m4_bmass)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m4.d <- Deriv(m4_bmass)
Term="month"
m4.dci <- confint(m4.d, term = Term)
m4.dsig <- signifD(pdat2$p3, d = m4.d[[Term]]$deriv,
                   +m4.dci[[Term]]$upper, m4.dci[[Term]]$lower)
plot.Deriv(m4.d, sizer=T, term=Term)
abline(v=4.1, type="l", lty=2)

#PP MALES####
PP_male_con=pp_plot%>%filter(treatment=="control", sex=="male")
PP_male_ex=pp_plot%>%filter(treatment=="exclosure", sex=="male")
PP_female_con=pp_plot%>%filter(treatment=="control", sex=="female")
PP_female_ex=pp_plot%>%filter(treatment=="exclosure", sex=="female")

#control####
#no explicit biotic effects
m5=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PP_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m5) #R2=57.3%

#add effect of DM+PB biomass
m5_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=PP_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m5_bmass) #R2=57.3%

#add effect of DM+PB biomass AND intraspecific competition (PP competition)
m5_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=PP_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m5_bmass2) #R2=59.6%

#best fit plot####
plot(m5_bmass2, shade=T, pages=1)
want=seq(1, nrow(PP_male_con), length.out = 200)
pdat2=with(PP_male_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                    precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                    bmass_PB=bmass_PB[want],bmass_PP= bmass_PP[want]))
p3=predict(m5_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m5_bmass2)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m5.d <- Deriv(m5_bmass2)
Term="month"
m5.dci <- confint(m5.d, term = Term)
m5.dsig <- signifD(pdat2$p3, d = m5.d[[Term]]$deriv,
                   +m5.dci[[Term]]$upper, m5.dci[[Term]]$lower)
plot.Deriv(m5.d, sizer=T, term=Term)
abline(v=3.3, type="l", lty=2)

#exclosure####

#no explicit biotic effects
m6=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m6) #R2=59.2%

#add effect of DM+PB biomass
m6_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m6_bmass) #R2=59.2%

#add effect of DM+PB biomass AND intraspecific competition (PB competition)
m6_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

summary(m6_bmass2) #R2=59.6%

#best fit plot####
plot(m6_bmass, shade=T, pages=1)
want=seq(1, nrow(PB_male_ex), length.out = 200)
pdat2=with(PB_male_ex, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                    precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                    bmass_PB=bmass_PB[want], bmass_PP=bmass_PP[want]))
p3=predict(m6_bmass, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m6_bmass)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m6.d <- Deriv(m6_bmass)
Term="month"
m6.dci <- confint(m6.d, term = Term)
m6.dsig <- signifD(pdat2$p3, d = m6.d[[Term]]$deriv,
                   +m6.dci[[Term]]$upper, m6.dci[[Term]]$lower)
plot.Deriv(m6.d, sizer=T, term=Term)
abline(v=4, type="l", lty=2)

#PP FEMALES####
#control####
#no explicit biotic effects
m7=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PP_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m7) #R2=46.6%

#add effect of DM+PB biomass
m7_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=PP_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m7_bmass) #R2=49%

#add effect of DM+PB biomass AND intraspecific competition (PB competition)
m7_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=PP_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m7_bmass2) #R2=52.3%

#best fit plot####
plot(m7_bmass2, shade=T, pages=1)
want=seq(1, nrow(PP_female_con), length.out = 200)
pdat2=with(PP_female_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                   precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                   bmass_PB=bmass_PB[want],bmass_PP= bmass_PP[want]))
p3=predict(m7_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m7_bmass2)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m7.d <- Deriv(m7_bmass2)
Term="month"
m7.dci <- confint(m7.d, term = Term)
m7.dsig <- signifD(pdat2$p3, d = m7.d[[Term]]$deriv,
                   +m7.dci[[Term]]$upper, m7.dci[[Term]]$lower)
plot.Deriv(m7.d, sizer=T, term=Term)
abline(v=5.1, type="l", lty=2)

#exclosure####
#no explicit biotic effects
m8=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m8) #R2=33.6%

#add effect of DM+PB biomass
m8_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m8_bmass) #R2=39.8%

#add effect of DM+PB biomass AND intraspecific competition (PB competition)
m8_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

summary(m8_bmass2) #R2=40.5%

#best fit plot####
plot(m8_bmass, shade=T, pages=1)
want=seq(1, nrow(PP_female_ex), length.out = 200)
pdat2=with(PP_female_ex, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                     precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                     bmass_PB=bmass_PB[want],bmass_PP= bmass_PP[want]))
p3=predict(m8_bmass, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m8_bmass)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m8.d <- Deriv(m8_bmass)
Term="month"
m8.dci <- confint(m8.d, term = Term)
m8.dsig <- signifD(pdat2$p3, d = m8.d[[Term]]$deriv,
                   +m8.dci[[Term]]$upper, m8.dci[[Term]]$lower)
plot.Deriv(m8.d, sizer=T, term=Term)
abline(v=5.4, type="l", lty=2)

#DM MALES####
DM_male_con=dm_plot%>%filter(treatment=="control", sex=="male")
DM_male_ex=dm_plot%>%filter(treatment=="exclosure", sex=="male")
DM_female_con=dm_plot%>%filter(treatment=="control", sex=="female")
DM_female_ex=dm_plot%>%filter(treatment=="exclosure", sex=="female")

#control####
#no explicit biotic effects
m9=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=DM_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m9) #R2=13.1%

#add effect of DM+PB biomass
m9_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=DM_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m9_bmass) #R2=14.4%

#add effect of DM+PB biomass AND intraspecific competition (PP competition)
m9_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=DM_male_con, method = 'REML', weights = abundance, family = binomial)

summary(m9_bmass2) #R2=17%

#best fit plot####
plot(m9_bmass2, shade=T, pages=1)
want=seq(1, nrow(DM_male_con), length.out = 200)
pdat2=with(DM_male_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                   precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                   bmass_PB=bmass_PB[want],bmass_PP= bmass_PP[want]))
p3=predict(m9_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m9_bmass2)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m9.d <- Deriv(m9_bmass2)
Term="month"
m9.dci <- confint(m9.d, term = Term)
m9.dsig <- signifD(pdat2$p3, d = m9.d[[Term]]$deriv,
                   +m9.dci[[Term]]$upper, m9.dci[[Term]]$lower)
plot.Deriv(m9.d, sizer=T, term=Term)
abline(v=3.2, type="l", lty=2)


#DM FEMALES####
#control####
#no explicit biotic effects
m10=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)
             , data=DM_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m10) #R2=27.2%

#add effect of DM+PB biomass
m10_bmass=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                     s(bmass_PB),
                   data=DM_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m10_bmass) #R2=27.9%

#add effect of DM+PB biomass AND intraspecific competition (PB competition)
m10_bmass2=mgcv::gam(proportion~s(month, bs="cc", k=12)+s(year)+s(ndvi)+s(precipitation)+s(bmass_DM)+
                      s(bmass_PB)+s(bmass_PP), data=DM_female_con, method = 'REML', weights = abundance, family = binomial)

summary(m10_bmass2) #R2=33.4%

#best fit plot####
plot(m10_bmass2, shade=T, pages=1)
want=seq(1, nrow(DM_female_con), length.out = 200)
pdat2=with(DM_female_con, data.frame(year=year[want], month=month[want], ndvi=ndvi[want], 
                                   precipitation=precipitation[want], bmass_DM=bmass_DM[want],
                                   bmass_PB=bmass_PB[want],bmass_PP= bmass_PP[want]))
p3=predict(m10_bmass2, newdata=pdat2, type="terms", se.fit = TRUE)
pdat2=transform(pdat2, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(m10_bmass2)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat2=transform(pdat2, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m10.d <- Deriv(m10_bmass2)
Term="month"
m10.dci <- confint(m10.d, term = Term)
m10.dsig <- signifD(pdat2$p3, d = m10.d[[Term]]$deriv,
                   +m10.dci[[Term]]$upper, m10.dci[[Term]]$lower)
plot.Deriv(m10.d, sizer=T, term=Term)
abline(v=4.3, type="l", lty=2)
