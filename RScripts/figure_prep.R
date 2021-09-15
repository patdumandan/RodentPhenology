#code for plotting#######
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

# plots##############

# PB abiotic only#####
pbf1=dwplot(list(pbf_con1, pbf_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Exclosure","Control"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm1=dwplot(list(pbm_con1, pbm_ex1),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Exclosure","Control"))+
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
                   labels = c("Exclosure","Control"))+
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
                   labels = c("Exclosure","Control"))+
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
                   labels = c("Exclosure","Control"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm3=dwplot(list(pbm_con3, pbm_ex3),vars_order = c("ndvis","ndvis_lag", "temps_mean","temps_lag_mean", 
                                                   "ppts_warm","ppts_lag_warm", 
                                                   "ppts_cool",  "ppts_lag_cool", "pbs", "pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",ndvis_lag="NDVI (lag 1)",
                       temps_mean="mean temperature (lag 0)", temps_lag_mean="mean temperature(lag 1)",
                       ppts_warm= "warm precipitation (lag 0)", ppts_lag_warm=" warm precipitation (lag 1)",
                       ppts_cool= "cool precipitation (lag 0)", ppts_lag_cool= "cool precipitation (lag 1)",
                       pbs="PB biomass", pps= "PP biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_grey(start = .3,end = .7,name = "Treatment",
                   labels = c("Exclosure","Control"))+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pb3=ggarrange(pbf3,pbm3, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))

annotate_figure(pb3, top=text_grob("abiotic + intra- + interspecific competition", face="italic", size=14, hjust=1))

