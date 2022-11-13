#this is script derived from G. Simpson's blog on idenifying periods of change

## download the derivatives gist
tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)

#load packages####

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dotwhisker)
library(ggpubr)
library(lattice)
library(mgcv)

#load data####
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

#models####
pbm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_con, method = 'REML', weights = abundance, family = binomial)
pbm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

pbf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_con, method = 'REML', weights = abundance, family = binomial)
pbf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

ppm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_con, method = 'REML', weights = abundance, family = binomial)
ppm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

ppf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_con, method = 'REML', weights = abundance, family = binomial)
ppf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

dipof_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipo_f, method = 'REML', weights = abundance, family = binomial)
dipom_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipo_m, method = 'REML', weights = abundance, family = binomial)

#DIPOS female####
want=seq(1, nrow(dipo_f), length.out = 200)
pdat=with(dipo_f, data.frame(year=year[want], month=month[want], 
                                    ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                    ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                    ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                    temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                    pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(dipof_con3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit, se3=p3$se.fit) #p2=fit, se2=std.error

df.res=df.residual(dipof_con3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(dipof_con3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

plot(dipo_f$proportion~dipo_f$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female control")

#male####
want=seq(1, nrow(dipo_m), length.out = 200)
pdat=with(dipo_m, data.frame(year=year[want], month=month[want], 
                             ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                             ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                             ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                             temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                             pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(dipom_con3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit, se3=p3$se.fit) #p2=fit, se2=std.error

df.res=df.residual(dipom_con3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(dipom_con3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#PB female####
#control/exclosure

want=seq(1, nrow(PB_female_ex), length.out = 200)
pdat=with(PB_female_ex, data.frame(year=year[want], month=month[want], 
                                    ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                    ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                    ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                    temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                    pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(pbf_ex3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit, se3=p3$se.fit) #p2=fit, se2=std.error

df.res=df.residual(pbf_ex3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(pbf_ex3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

plot(PB_female_con$proportion~PB_female_con$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female control")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)
rect(xleft=11, ybottom=0, ytop=1, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=0, ytop=1, xright=1.5, col = rgb(0,0,0.5, 1/4))#control
rect(xleft=5, ybottom=0, ytop=1, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

plot(PB_female_ex$proportion~PB_female_ex$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female exclosure")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)
rect(xleft=8, ybottom=0, ytop=1, xright=12, col = rgb(0,0,0.5, 1/4)) #exclosure incr
rect(xleft=1, ybottom=0, ytop=1, xright=1.5, col = rgb(0,0,0.5, 1/4)) #exclosure incr
rect(xleft=5, ybottom=0, ytop=1, xright=8, col = rgb(0.5,0,0, 1/4))#exclosure decr


#PB male####
plot(PB_male_con$proportion~PB_male_con$month, xaxt="n", xlab="month", ylab="P(breeding)")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)


want=seq(1, nrow(PB_male_ex), length.out = 200)
pdat=with(PB_male_ex, data.frame(year=year[want], month=month[want], 
                                 ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                 ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                 ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                 temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                 pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(pbm_ex3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(pbm_ex3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(pbm_ex3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
abline(v=3.1, type="l", lty=2)

#PP female####
#control exclosure

want=seq(1, nrow(PP_female_ex), length.out = length(PP_female_ex$X))
pdat=with(PP_female_ex, data.frame(year=year[want], month=month[want], 
                                    ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                    ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                    ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                    temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                    pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(ppf_ex3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(ppf_ex3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(ppf_ex3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

plot(PP_female_ex$proportion~PP_female_ex$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female exclosure")

want=seq(1, nrow(PP_female_con), length.out = 200)
pdat=with(PP_female_con, data.frame(year=year[want], month=month[want], 
                                   ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                   ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                   ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                   temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                   pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(ppf_con3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(ppf_con3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(ppf_con3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

plot(PP_female_ex$proportion~PP_female_ex$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female exclosure")




