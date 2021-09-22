tot_dat=read.csv("./reconfigured_data/raw_cleaned/reprod_propn_male.csv")
require(dplyr)
PB_dat_M=tot_dat%>%filter(species=="PB", !(treatment=="spectabs"))
PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PB_dat_M$mon_cos= cos(2*pi*(PB_dat_M$month/12))
PB_dat_M$mon_sin= sin(2*pi*(PB_dat_M$month/12))
PB_dat_M$spcode=as.integer(PB_dat_M$species)

ggplot(PB_dat_M, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB males")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#for seasonal gam
PB_gam_con=PB_dat_M%>%
  filter(treatment=="control")%>%
  mutate(date=as.POSIXct(paste(month, "15", year, sep="-"), format="%m-%d-%Y"), Time=as.numeric(date)/1000,
         monthtxt=lubridate::month(date, label=T))

PB_gam_ex=PB_dat_M%>%
  filter(treatment=="exclosure")%>%
  mutate(date=as.POSIXct(paste(month, "15", year, sep="-"), format="%m-%d-%Y"), Time=as.numeric(date)/1000,
         monthtxt=lubridate::month(date, label=T))
  

require(mgcv)

m1=gamm(proportion~s(month, bs="cc", k=12), data=PB_gam_con, correlation = corARMA(form= ~1|year, p=2))
m2=gamm(proportion~s(month, bs="cc", k=12), data=PB_gam_ex, correlation = corARMA(form= ~1|year, p=2))
m1
str(PB_gam_con)
plot(m1$gam)
points(PB_gam_ex$proportion)
plot(PB_gam_con$proportion~PB_gam_con$month)
plot(PB_gam_ex$proportion~PB_gam_ex$month)
length(unique(PB_gam_con$year))
m1pred=predict(m1$gam)
plot(m1$gam)

#CALCULATE DERIVATIVES

tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)
ls()

#control
want=seq(1, nrow(PB_female_con), length.out = 200)
pdat=with(PB_female_con, data.frame(year=year[want], month=month[want], 
                                    ndvis=ndvis[want],ndvis_lag=ndvis_lag[want],
                                    ppts_warm=ppts_warm[want], ppts_lag_warm=ppts_lag_warm[want],
                                    ppts_cool=ppts_cool[want], ppts_lag_cool=ppts_lag_cool[want],
                                    temps_mean=temps_mean[want], temps_lag_mean=temps_lag_mean[want],
                                    pps=pps[want], pbs=pbs[want], dipos=dipos[want]))
p3=predict(pbf_con3, newdata=pdat, type="terms", se.fit = TRUE)
pdat=transform(pdat, p3=p3$fit[,1], se3=p3$se.fit[,1]) #p2=fit, se2=std.error

df.res=df.residual(pbf_con3)

crit.t=qt(0.025, df.res, lower.tail = F)
pdat=transform(pdat, upper=p3+(crit.t*se3), lower=p3-(crit.t*se3))

m1.d <- Deriv(pbf_con3)
Term="month"
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p3, d = m1.d[[Term]]$deriv,
                   +m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)
plot.Deriv(m1.d, sizer=T, term=Term, xaxt="n")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
abline(v=3.1, type="l", lty=2)

ggplot(PB_dat_M, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB males")+
  theme_classic()
