#code for plotting#######
library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(dotwhisker)
library(ggpubr)

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

#pb models####
pbm_con1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PB_male_con, method = 'REML', weights = abundance, family = binomial)
pbm_ex1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

pbf_con1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PB_female_con, method = 'REML', weights = abundance, family = binomial)
pbf_ex1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

pbm_con2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs, data=PB_male_con, method = 'REML', weights = abundance, family = binomial)
pbm_ex2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs, data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

pbf_con2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs, data=PB_female_con, method = 'REML', weights = abundance, family = binomial)
pbf_ex2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs, data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)

pbm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_con, method = 'REML', weights = abundance, family = binomial)
pbm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_male_ex, method = 'REML', weights = abundance, family = binomial)

pbf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_con, method = 'REML', weights = abundance, family = binomial)
pbf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PB_female_ex, method = 'REML', weights = abundance, family = binomial)
#phenology models####

pbplot=ggplot(pb_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Exclosure","Control"))+theme_pubr(base_size = 10, x.text.angle =90)
  

ppplot=ggplot(pp_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
   facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Exclosure","Control"))+theme_pubr(base_size = 10, x.text.angle=90)

ggarrange(pbplot,ppplot, nrow=1, common.legend = T, legend="right",
              labels=c("A", "B"))

# drivers plots##############

# PB abiotic only#####
pbf1=dwplot(list(pbf_con1, pbf_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm1=dwplot(list(pbm_con1, pbm_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


pb1=ggarrange(pbf1,pbm1, nrow=1, common.legend = T, legend="right",
              labels=c("female", "male"))

annotate_figure(pb1, top=text_grob("abiotic only", face="italic", size=14, hjust=3))

#PB abiotic + intraspecific####
pbf2=dwplot(list(pbf_con2, pbf_ex2),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs= "PB biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm2=dwplot(list(pbm_con2, pbm_ex2),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="PB biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pb2=ggarrange(pbf2,pbm2, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))

annotate_figure(pb2, top=text_grob("abiotic + intraspecific competition", face="italic", size=14, hjust=1))

#PB abiotic + intra + interspecific competition####
pbf3=dwplot(list(pbf_con3, pbf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs", "pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs= "PB biomass", pps= "PP biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle("A")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbf3=pbf3%>% annotate_figure(top=text_grob("female", face = "bold"))

pbm3=dwplot(list(pbm_con3, pbm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs", "pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="PB biomass", pps= "PP biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle("B")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm3=pbm3%>% annotate_figure(top=text_grob("male", face = "bold"))

pb3=ggarrange(pbf3,pbm3, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))

annotate_figure(pb3, top=text_grob("Bailey's pocket mouse", face="italic", size=14, hjust=1))

#pp models####
ppm_con1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PP_male_con, method = 'REML', weights = abundance, family = binomial)
ppm_ex1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

ppf_con1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PP_female_con, method = 'REML', weights = abundance, family = binomial)
ppf_ex1=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

ppm_con2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pps, data=PP_male_con, method = 'REML', weights = abundance, family = binomial)
ppm_ex2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pps, data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

ppf_con2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pps, data=PP_female_con, method = 'REML', weights = abundance, family = binomial)
ppf_ex2=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pps, data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)

ppm_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_con, method = 'REML', weights = abundance, family = binomial)
ppm_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_male_ex, method = 'REML', weights = abundance, family = binomial)

ppf_con3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_con, method = 'REML', weights = abundance, family = binomial)
ppf_ex3=mgcv::gam(proportion~s(month,bs="cc")+s(year)+ndvis+ndvis_lag+temps_mean+temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool+ pbs+pps+dipos, data=PP_female_ex, method = 'REML', weights = abundance, family = binomial)
#PP abiotic####

ppf1=dwplot(list(ppf_con1, ppf_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm1=dwplot(list(ppm_con1, ppm_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


pp1=ggarrange(ppf1,ppm1, nrow=1, common.legend = T, legend="right",
              labels=c("female", "male"))

annotate_figure(pp1, top=text_grob("abiotic only", face="italic", size=14, hjust=3))

#PP abiotic + intraspecific####
ppf2=dwplot(list(ppf_con2, ppf_ex2),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pps"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pps= "PP biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm2=dwplot(list(ppm_con2, ppm_ex2),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pps"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pps="PP biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pp2=ggarrange(ppf2,ppm2, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))

annotate_figure(pp2, top=text_grob("abiotic + intraspecific competition", face="italic", size=14, hjust=1))

#PP abiotic + intra +interspecific####
ppf3=dwplot(list(ppf_con3, ppf_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pps","pbs", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                      pps= "PP biomass",pbs= "PB biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle ("C")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm3=dwplot(list(ppm_con3, ppm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pps","pbs", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pps= "PP biomass", pbs="PB biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Control","Exclosure"))+ggtitle("D")+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pp3=ggarrange(ppf3,ppm3, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))

annotate_figure(pp3, top=text_grob("desert pocket mouse", face="italic", size=14, hjust=1))

#phenology plots####

