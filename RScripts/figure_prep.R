#code for plotting#######
library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dotwhisker)
library(ggpubr)
library(lattice)

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

#GAMs####

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

#FIGURE 1####

#DIPO female control####
plot.gam(dipof_con3, select=1, xaxt="n", ylab="Effect")
title("female")
mtext("kangaroo rat", at=1, font=2)
mtext("A", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=11.5, ybottom=-5, ytop=5, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#DIPO male control####
plot(dipo_m$proportion~dipo_m$month, xaxt="n", pch=20, ylab="proportion of breeders", xlab="month")
title("male")
#plot.gam(dipom_con3, select=1, xaxt="n", ylab="Effect", main="male")
axis(1, at=1:12, labels=F)
mtext("B", at=12, col="red", font=2)

#PB female control####

plot.gam(pbf_con3, select=1, xaxt="n", ylab="Effect")
mtext("Bailey's pocket mouse (control)", at=3, font=2)
mtext("C", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PB male control####
plot(PB_male_con$proportion~PB_male_con$month, xaxt="n", pch=20, ylab="proportion of breeders", xlab="month")
mtext("D", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)

#PB female exclosure####
plot.gam(pbf_ex3, select=1, xaxt="n", ylab="Effect")
mtext("Bailey's pocket mouse (exclosure)", at=3, font=2)
mtext("E", at=12, col="red", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=2, ybottom=-5, ytop=5, xright=4, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PB male exclosure####
plot.gam(pbm_ex3, select=1, xaxt="n", ylab="Effect")
mtext("F", at=12, col="red", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=7, ybottom=-5, ytop=5, xright=7.5, col = rgb(0.5,0,0, 1/4))#control decr

#FIGURE 2####
#DIPO female control####
plot.gam(dipof_con3, select=1, xaxt="n", ylab="Effect")
title("female")
mtext("kangaroo rat", at=1, font=2)
mtext("A", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=11.5, ybottom=-5, ytop=5, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#DIPO male control####
plot(dipo_m$proportion~dipo_m$month, xaxt="n", pch=20, ylab="proportion of breeders", xlab="month")
title("male")
#plot.gam(dipom_con3, select=1, xaxt="n", ylab="Effect", main="male")
axis(1, at=1:12, labels=F)
mtext("B", at=12, col="red", font=2)

#PP female control####

plot.gam(ppf_con3, select=1, xaxt="n", ylab="Effect")
mtext("desert pocket mouse (control)", at=3, font=2)
mtext("C", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=11, ybottom=-5, ytop=10, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=-5, ytop=10, xright=1.5, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=5, ybottom=-5, ytop=10, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PP male control####
plot.gam(ppm_con3, select=1, xaxt="n", ylab="Effect")
mtext("D", at=12, col="red", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=1, ybottom=-5, ytop=10, xright=4, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=10, xright=9, col = rgb(0.5,0,0, 1/4))#control decr

#PP female exclosure####

plot.gam(ppf_ex3, select=1, xaxt="n", ylab="Effect")
mtext("desert pocket mouse (exclosure)", at=3, font=2)
mtext("E", at=12, col="red", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=10, ybottom=-5, ytop=5, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=5, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr


#PP male exclosure####
plot.gam(ppm_ex3, select=1, xaxt="n", ylab="Effect")
mtext("F", at=12, col="red", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=1, ybottom=-5, ytop=5, xright=2, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=7, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#FIGURE 3####

#PB####
pbd_f_con=dwplot(list(dipof_con3,pbf_con3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                               "ppts_warm","ppts_lag_warm", 
                                                                               "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green"),
                     labels = c("kangaroo rats","Bailey's pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green"),
                    labels = c("kangaroo rats","Bailey's pocket mice"))+
  ggtitle("female")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbd_m_con=dwplot(list(dipom_con3,pbm_con3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                               "ppts_warm","ppts_lag_warm", 
                                                                               "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green"),
                     labels = c("kangaroo rats","Bailey's pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green"),
                    labels = c("kangaroo rats","Bailey's pocket mice"))+
  ggtitle("male")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pb1=ggarrange(pbd_f_con, pbd_m_con+ theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
annotate_figure(pb1, top=text_grob("control", face="bold"))

#exclosure####

#female####
pbpp_f_ex=dwplot(list(pbf_ex3,ppf_ex3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                           "ppts_warm","ppts_lag_warm", 
                                                                           "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("green", "orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("female")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

#male####
pbpp_m_ex=dwplot(list(pbm_ex3,ppm_ex3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                           "ppts_warm","ppts_lag_warm", 
                                                                           "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("green", "orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("male")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbpp=ggarrange(pbpp_f_ex, pbpp_m_ex+theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
               font.label = list(size=12))

annotate_figure(pbpp, top=text_grob("exclosure", face="bold"))

#FIGURE 4####

#PP####
ppd_f_con=dwplot(list(dipof_con3,ppf_con3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                               "ppts_warm","ppts_lag_warm", 
                                                                               "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "orange"),
                     labels = c("kangaroo rats","desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "orange"),
                    labels = c("kangaroo rats","desert pocket mice"))+
  ggtitle("female")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppd_m_con=dwplot(list(dipom_con3,ppm_con3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                               "ppts_warm","ppts_lag_warm", 
                                                                               "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "orange"),
                     labels = c("kangaroo rats","desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "orange"),
                    labels = c("kangaroo rats","desert pocket mice"))+
  ggtitle("male")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pp1=ggarrange(ppd_f_con, ppd_m_con+ theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
annotate_figure(pp1, top=text_grob("control", face="bold"))

#exclosure####

#female####
pbpp_f_ex=dwplot(list(pbf_ex3,ppf_ex3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                           "ppts_warm","ppts_lag_warm", 
                                                                           "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("green", "orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("female")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

#male####
pbpp_m_ex=dwplot(list(pbm_ex3,ppm_ex3),style="distribution",vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                                           "ppts_warm","ppts_lag_warm", 
                                                                           "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("green", "orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("male")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbpp=ggarrange(pbpp_f_ex, pbpp_m_ex+theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
               font.label = list(size=12))

annotate_figure(pbpp, top=text_grob("exclosure", face="bold"))



#APPENDIX FIGURE 3####

#PB####
pbf3=dwplot(list(ppf_con3, ppf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs= "Bailey's pocket mouse biomass", pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(values=c( "darkblue", "brown"),name = "Treatment",
                     labels = c("Control","Exclosure"))+ggtitle ("C")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm3=dwplot(list(ppm_con3, ppm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(values=c( "darkblue", "brown"),name = "Treatment",
                     labels = c("Control","Exclosure"))+ggtitle("D")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pb3=ggarrange(pbf3,pbm3+theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))

#PP ####
ppf3=dwplot(list(ppf_con3, ppf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs= "Bailey's pocket mouse biomass", pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(values=c( "darkblue", "brown"),name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle ("C")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm3=dwplot(list(ppm_con3, ppm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(values=c( "darkblue", "brown"),name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle("D")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pp3=ggarrange(ppf3,ppm3+theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))

#annotate_figure(pp3, top=text_grob("desert pocket mouse", face="italic", size=14, hjust=1))

#DIPO####
dipo3=dwplot(list(dipof_con3, dipom_con3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                         "ppts_warm","ppts_lag_warm", 
                                                         "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_discrete(name = "Sex",
                   labels = c("Female","Male"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)
