#PB female####
plot(PP_female_con$proportion~PP_female_con$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female control")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)
rect(xleft=11, ybottom=0, ytop=1, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=0, ytop=1, xright=1.5, col = rgb(0,0,0.5, 1/4))#control
rect(xleft=5, ybottom=0, ytop=1, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

plot(PP_female_ex$proportion~PP_female_ex$month, xaxt="n", xlab="month", ylab="P(breeding)", pch=16, main="PP female exclosure")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)
rect(xleft=8, ybottom=0, ytop=1, xright=12, col = rgb(0,0,0.5, 1/4)) #exclosure incr
rect(xleft=1, ybottom=0, ytop=1, xright=1.5, col = rgb(0,0,0.5, 1/4)) #exclosure incr
rect(xleft=5, ybottom=0, ytop=1, xright=8, col = rgb(0.5,0,0, 1/4))#exclosure decr


#PB male####
plot(PP_male_con$proportion~PP_male_con$month, xaxt="n", xlab="month", ylab="P(breeding)")
axis(1, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"), at=c(1:12), srt=90)

rect(xleft=1, ybottom=0, ytop=1, xright=4, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=0, ytop=1, xright=2, col = rgb(0,0.5,0, 1/4)) #exclosure incr

rect(xleft=6, ybottom=0, ytop=1, xright=9, col = rgb(0.5,0,0, 1/4))#control decr
rect(xleft=7, ybottom=0, ytop=1, xright=8, col = rgb(0,0.5,0.5, 1/4))#exclosure decr


want=seq(1, nrow(PP_male_ex), length.out = 200)
pdat=with(PP_male_ex, data.frame(year=year[want], month=month[want], 
                                 ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                 ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                 ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                 temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                 pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(ppm_ex3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(ppm_ex3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(ppm_ex3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
abline(v=3.1, type="l", lty=2)

#PB MALE

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
