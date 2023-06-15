#code for plotting##
library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dotwhisker)
library(ggpubr)
library(lattice)

#FIGURE 1####

#PB models####

pbf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbf_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n",main="female", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

pbm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbm_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", main="male",ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#PP models####

ppf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppf_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

ppm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppm_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):k-rat inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#FIGURE 2####

#DIPO female control####
plot.gam(dipof_con3, select=1, xaxt="n")
title("female")
#text(x=11, y=1, "(a)", col="black", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=11.5, ybottom=-5, ytop=5, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#DIPO male control####
plot(dipo_m$proportion~dipo_m$month, xaxt="n", pch=20, ylab="proportion of breeders")
title("male")
#plot.gam(dipom_con3, select=1, xaxt="n", ylab="Effect", main="male")
#axis(1, at=1:12, labels=F)
#text(x=11, y=0.9, "(b)", col="black", font=2)

#PB female control####

plot.gam(pbf_con3, select=1, xaxt="n")
#mtext("Bailey's pocket mouse (control)", at=3, font=2)
#text(x=11, y=3.1, "(c)", col="black", font=2)
axis(1, at=1:12, labels=F)
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PB male control####
plot(PB_male_con$proportion~PB_male_con$month, xaxt="n", pch=20, ylab="proportion of breeders", xlab="month")
#text(x=11, y=0.9, "(d)", col="black", font=2)
axis(1, at=1:12, labels=F)

#PB female exclosure####
plot.gam(pbf_ex3, select=1, xaxt="n")
#mtext("Bailey's pocket mouse (exclosure)", at=3, font=2)
#text(x=11, y=2.7, "(e)", col="black", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=2, ybottom=-5, ytop=5, xright=4, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PB male exclosure####
plot.gam(pbm_ex3, select=1, xaxt="n")
#text(x=11, y=1.3, "(f)", col="black", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=1, ybottom=-5, ytop=5, xright=3, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=7, ybottom=-5, ytop=5, xright=7.5, col = rgb(0.5,0,0, 1/4))#control decr

#FIGURE 3####

#PB ####
pbf=dwplot(list(pbf_con3,pbf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                 "ppts_warm","ppts_lag_warm", 
                                                 "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  ggtitle("female")+#geom_label(label="A", x=3, y=11, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm=dwplot(list(pbm_con3,pbm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                 "ppts_warm","ppts_lag_warm", 
                                                 "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  ggtitle("male")+#geom_label(label="B", x=3, y=11, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


pb1=ggarrange(pbf, pbm + theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
pb1
annotate_figure(pb1, top=text_grob("Bailey's pocket mouse", face="bold"))

#PP####
ppf=dwplot(list(ppf_con3,ppf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                 "ppts_warm","ppts_lag_warm", 
                                                 "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","exclosure"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","exclosure"))+
  geom_label(label="C", x=3, y=11, label.size = NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm=dwplot(list(ppm_con3,ppm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                 "ppts_warm","ppts_lag_warm", 
                                                 "ppts_cool",  "ppts_lag_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","exclosure"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","exclosure"))+
  geom_label(label="D", x=1.5, y=11, label.size = NA)+xlim(-4,4)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pp1=ggarrange(ppf, ppm + theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
pp1
annotate_figure(pp1, top=text_grob("desert pocket mouse", face="bold"))

#FIGURE 4####

#PP female control####

plot.gam(ppf_con3, select=1, xaxt="n")
#text(x=11, y=4, "(c)", col="black", font=2)
rect(xleft=11, ybottom=-5, ytop=10, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=1, ybottom=-5, ytop=10, xright=1.5, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=5, ybottom=-5, ytop=10, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#PP male control####
plot.gam(ppm_con3, select=1, xaxt="n")
#text(x=11, y=2.4, "(d)", col="black", font=2)
rect(xleft=1, ybottom=-5, ytop=10, xright=4, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=6, ybottom=-5, ytop=10, xright=9, col = rgb(0.5,0,0, 1/4))#control decr

#PP female exclosure####

plot.gam(ppf_ex3, select=1, xaxt="n")
#text(x=11, y=3.4, "(e)", col="black", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=10, ybottom=-5, ytop=5, xright=12, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=5, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr


#PP male exclosure####
plot.gam(ppm_ex3, select=1, xaxt="n")
#text(x=11, y=1.4, "(f)", col="black", font=2)
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))
rect(xleft=1, ybottom=-5, ytop=5, xright=2, col = rgb(0,0,0.5, 1/4)) #control
rect(xleft=7, ybottom=-5, ytop=5, xright=8, col = rgb(0.5,0,0, 1/4))#control decr

#FIGURE 5####
#control####
plot_con1=dwplot(list(dipof_con3,pbf_con3, pbf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("female")+geom_label(label="A", x=4, y=11, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


plot_con2=dwplot(list(dipom_con3,pbm_con3, ppm_con3),style="distribution",vars_order = c("ndvis", "temps_mean",
                                                                                         "ppts_warm", 
                                                                                         "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("male")+ geom_label(label="B", x=2, y=13, label.size=NA)+
  geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con=ggarrange(plot_con1, plot_con2+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                   font.label = list(size=12))
annotate_figure(plot_con, top=text_grob("control", face="bold"))
#exclosure####

plot_ex1=dwplot(list(pbf_ex3, ppf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                           "ppts_warm",
                                                                           "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c( "green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c( "green","orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("female")+geom_label(label="C", x=4, y=9.5, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


plot_ex2=dwplot(list(pbm_ex3, ppm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean",
                                                                             "ppts_warm", 
                                                                             "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c( "green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c( "green","orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("male")+ geom_label(label="D", x=2, y=9, label.size=NA)+
  geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_ex=ggarrange(plot_ex1, plot_ex2+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                  font.label = list(size=12))
annotate_figure(plot_ex, top=text_grob("k-rat inaccessible", face="bold"))

