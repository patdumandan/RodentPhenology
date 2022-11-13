all_dat=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/reproductive_full_data.csv")

DM_all=all_dat%>%filter(species=="DM")
DO_all=all_dat%>%filter(species=="DO")
DS_all=all_dat%>%filter(species=="DS")
Dipos_all=all_dat%>%filter(species%in%c("DM", "DO", "DS"))

dipos_female=Dipos_all%>%filter(sex=="female")  
dipos_male=Dipos_all%>%filter(sex=="male")  

dipo_f3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipos_female, method = 'REML', weights = abundance, family = binomial)

dipo_m3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=dipos_male, method = 'REML', weights = abundance, family = binomial)

want=seq(1, nrow(dipos_male), length.out = 200)
pdat=with(dipos_male, data.frame(year=year[want], month=month[want], 
                                    ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                    ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                    ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                    temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                    pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(dipo_m3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(dipo_m3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(dipo_m3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
annotate("text", label="Dipo female")
