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

#FIGURE 2####

#PB ####
pbf=dwplot(list(pbf_con3,pbf_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm",
                                                 "ppts_cool",  "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)",
                       ppts_cool= "cool precipitation (lag 0)",
                       pbs="C. baileyi biomass",pps= "C. penicillatus biomass", dipos= "Dipodomys spp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("Dipodomys accessible","Dipodomys inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("Dipodomys accessible","Dipodomys inaccessible"))+
  ggtitle("female")+#geom_label(label="A", x=3, y=11, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm=dwplot(list(pbm_con3,pbm_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm",
                                                 "ppts_cool",  "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)",
                       ppts_cool= "cool precipitation (lag 0)",
                       pbs="C. baileyi biomass",pps= "C. penicillatus biomass", dipos= "Dipodomys spp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("Dipodomys accessible","Dipodomys inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("Dipodomys accessible","Dipodomys inaccessible"))+
  ggtitle("male")+#geom_label(label="A", x=3, y=11, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pb1=ggarrange(pbf, pbm + theme(axis.text.y = element_blank()),widths = c(0.35, 0.25), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
pb1

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

#FIGURE 3####

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

#FIGURE 4####
#PB plots####
plot_con11=dwplot(list(dipof_con3,pbf_con3, pbf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("skyblue", "green", "orange"),
                     labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("skyblue", "green","orange"),
                    labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  ggtitle("female")+#geom_label(label="A", x=4, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


plot_con12=dwplot(list(dipom_con3,pbm_con3, pbm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("skyblue", "green", "orange"),
                     labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("skyblue", "green","orange"),
                    labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  ggtitle("male")+#geom_label(label="A", x=4, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con=ggarrange(plot_con11, plot_con12+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                   font.label = list(size=12))
annotate_figure(plot_con, top=text_grob("k-rat inaccessible", face="bold"))

#PP plots####
plot_con21=dwplot(list(dipof_con3,ppf_con3, ppf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("skyblue", "green", "orange"),
                     labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("skyblue", "green","orange"),
                    labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  #ggtitle("female")+#geom_label(label="A", x=4, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


plot_con22=dwplot(list(dipom_con3,ppm_con3, ppm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("skyblue", "green", "orange"),
                     labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("skyblue", "green","orange"),
                    labels = c("kangaroo rats","Dipodomys accessible", "Dipodomys inaccessible"))+
  #ggtitle("male")+#geom_label(label="A", x=4, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con2=ggarrange(plot_con21, plot_con22+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                    font.label = list(size=12))
