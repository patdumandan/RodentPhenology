#APPENDIX####

#FIG 1####
use_default_data_path("D:\\Dropbox (UFL)\\PhD-stuff\\ReproPhenology")

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
  level="Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "crosstab",
  time = "all",
  output = "abundance",
  na_drop = FALSE,
  zero_drop = FALSE,
  min_traps = 1,
  min_plots = 1,
  effort = TRUE,
  download_if_missing = TRUE,
  quiet = FALSE
)%>%mutate(year=year(censusdate), month=month(censusdate))

pbcont_dat=rodent_data%>%
  select(year,treatment, PB)%>%
  filter(treatment=="control")%>%
  group_by(year)%>%
  summarize(total=sum(PB))

pbex_dat=rodent_data%>%
  select(year,treatment, PB)%>%
  filter(treatment=="exclosure")%>%
  group_by(year)%>%
  summarize(total=sum(PB))

pb1=ggplot(data=pbcont_dat, aes(year, total)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("k-rat accessible")+annotate("rect", alpha=.2, fill='red',xmin=1988, xmax=2014, ymin=0, ymax=480)

pb2=ggplot(data=pbex_dat, aes(year, total)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("k-rat inaccessible")+annotate("rect", alpha=.2, fill='red',xmin=1988, xmax=2014, ymin=0, ymax=480)

p2=ggarrange(pb1, pb2, common.legend = TRUE, ncol=2, nrow=1)
annotate_figure(p2, top = text_grob("C.baileyi", 
                                    face = "bold", size = 14))

#PP####

use_default_data_path("D:\\Dropbox (UFL)\\PhD-stuff\\ReproPhenology")

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
  level="Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "crosstab",
  time = "all",
  output = "abundance",
  na_drop = FALSE,
  zero_drop = FALSE,
  min_traps = 1,
  min_plots = 1,
  effort = TRUE,
  download_if_missing = TRUE,
  quiet = FALSE
)%>%mutate(year=year(censusdate), month=month(censusdate))

ppcont_dat=rodent_data%>%
  select(year,treatment, PP)%>%
  filter(treatment=="control")%>%
  group_by(year)%>%
  summarize(total=sum(PP))

ppex_dat=rodent_data%>%
  select(year,treatment, PP)%>%
  filter(treatment=="exclosure")%>%
  group_by(year)%>%
  summarize(total=sum(PP))

pp1=ggplot(data=ppcont_dat, aes(year, total)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("k-rat accessible")+annotate("rect", alpha=.2, fill='red',xmin=1988, xmax=2014, ymin=0, ymax=350)

pp2=ggplot(data=ppex_dat, aes(year, total)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("k-rat inaccessible")+annotate("rect", alpha=.2, fill='red',xmin=1988, xmax=2014, ymin=0, ymax=350)

p1=ggarrange(pp1, pp2, common.legend = TRUE, ncol=2, nrow=1)
annotate_figure(p1, top = text_grob("C.penicillatus", 
                                    face = "bold", size = 14))
#FIG 2####

#PB models####

pbf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbf_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n",main="female", ylab="s(difference):Dipodomys inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

pbm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+s(year)+
                 otrt, data=pbm_plot, method = 'REML', weights = abundance, family = binomial)

plot(pbm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", main="male",ylab="s(difference):Dipodomys inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#PP models####

ppf3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppf_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppf3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):Dipodomys inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

ppm3=mgcv::gam(proportion~s(month,bs="cc")+s(month,bs="cc", by= otrt)+otrt, data=ppm_plot, method = 'REML', weights = abundance, family = binomial)

plot(ppm3, shade = TRUE, scale = 0, seWithMean = TRUE, select=2, xaxt="n", ylab="s(difference):Dipodomys inaccessible")
axis(1, at=1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))

#FIG 3####

#PP###
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


#exclosure###

#female###
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

#male###
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

annotate_figure(pbpp, top=text_grob("k-rat inaccessible", face="bold"))

#PB###
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

#PP ###
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

#DIPO###
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

#Fig 3 alts###
#1###
#PB###
pbf=dwplot(list(pbf_con3,pbf_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm",
                                                 "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)",
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  ggtitle("female")+#geom_label(label="A", x=3, y=7, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

pbm=dwplot(list(pbm_con3,pbm_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm", 
                                                 "ppts_cool",  "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  ggtitle("male")+#geom_label(label="B", x=3, y=7, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


pb1=ggarrange(pbf, pbm + theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
pb1
annotate_figure(pb1, top=text_grob("Bailey's pocket mouse", face="bold"))

ppf=dwplot(list(ppf_con3,ppf_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm",
                                                 "ppts_cool", "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)",
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  #ggtitle("female")+#geom_label(label="A", x=3, y=7, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

ppm=dwplot(list(ppm_con3,pbm_ex3),vars_order = c("ndvis", "temps_mean", 
                                                 "ppts_warm", 
                                                 "ppts_cool",  "pbs","pps", "dipos"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Treatment",
                     values=c("skyblue", "green"),
                     labels = c("k-rat accessible","k-rat inaccessible"))+
  scale_fill_manual(name = "Treatment",
                    values=c("skyblue", "green"),
                    labels = c("k-rat accessible","k-rat inaccessible"))+
  #ggtitle("male")+#geom_label(label="B", x=3, y=7, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


pp1=ggarrange(ppf, ppm + theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
              font.label = list(size=12))
pp1
annotate_figure(pp1, top=text_grob("desert pocket mouse", face="bold"))

#2###

plot_ex1=dwplot(list(pbf_ex3, ppf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                           "ppts_warm",
                                                                           "ppts_cool", "pbs","pps"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)",
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass"))+
  scale_color_manual(name = "Species",
                     values=c( "green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c( "green","orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("female")+geom_label(label="C", x=4, y=8.75, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)


plot_ex2=dwplot(list(pbm_ex3, ppm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean",
                                                                           "ppts_warm", 
                                                                           "ppts_cool", "pbs","pps"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)", 
                       pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass"))+
  scale_color_manual(name = "Species",
                     values=c( "green", "orange"),
                     labels = c("Bailey's pocket mice", "desert pocket mice"))+
  scale_fill_manual(name = "Species",
                    values=c( "green","orange"),
                    labels = c("Bailey's pocket mice", "desert pocket mice"))+
  ggtitle("male")+ geom_label(label="D", x=2, y=8, label.size=NA)+
  geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_ex=ggarrange(plot_ex1, plot_ex2+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                  font.label = list(size=12))
annotate_figure(plot_ex, top=text_grob("k-rat inaccessible", face="bold"))

#PB regimes####
#female####
pf=pp_plot%>%filter(sex=="female")

p1=pf%>%filter(year%in%c(1977:1994))
a1=ggplot(p1, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("1977-1994")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA1=a1+geom_text(x=11, y=1, col="black", label="A")

p2=pf%>%filter(year%in%c(1995:2010))
a2=ggplot(p2, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("1995-2010")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA2=a2+geom_text(x=11, y=1, col="black", label="B")

p3=pf%>%filter(year%in%c(2011:2014))
a3=ggplot(p3, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("2011-2014")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA3=a3+geom_text(x=11, y=0.57, col="black", label="C")

a11=ggarrange(AA1,AA2,AA3, nrow=1, common.legend = T, legend="bottom")
annotate_figure(a11, top=text_grob("female", face="bold"))

#male####
pm=pp_plot%>%filter(sex=="male")

p12=pm%>%filter(year%in%c(1977:1994))
a12=ggplot(p12, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("1977-1994")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA12=a12+geom_text(x=11, y=1, col="black", label="D")

p22=pm%>%filter(year%in%c(1995:2010))
a22=ggplot(p22, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("1995-2010")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA22=a22+geom_text(x=11, y=1, col="black", label="E")

p32=pm%>%filter(year%in%c(2011:2014))
a32=ggplot(p32, aes(x=month, y=proportion, col=treatment))+geom_point()+theme_classic()+ggtitle("2011-2014")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  geom_smooth(method='gam')
AA32=a32+geom_text(x=11, y=0.73, col="black", label="F")

a22=ggarrange(AA12,AA22,AA32, nrow=1, common.legend = T, legend="bottom")
annotate_figure(a22, top=text_grob("male", face="bold"))

#new fig4 and 5####

#abiotic factors####

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

plot_con12=dwplot(list(dipof_con3,pbm_con3, pbm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  ggtitle("male")+#geom_label(label="B", x=2, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con1=ggarrange(plot_con11, plot_con12+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                    font.label = list(size=12))
annotate_figure(plot_con1, top=text_grob("Bailey's pocket mouse", face="bold"))

#PP plots####

plot_con21=dwplot(list(dipof_con3,ppf_con3, ppf_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","k-rat accessible", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","k-rat accessible", "k-rat inaccessible"))+
  #ggtitle("female")+#geom_label(label="C", x=4, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con22=dwplot(list(dipof_con3,ppm_con3, ppm_ex3),style="distribution",vars_order = c("ndvis", "temps_mean", 
                                                                                         "ppts_warm",
                                                                                         "ppts_cool"))%>%
  relabel_predictors(c(ndvis="NDVI (lag 0)",
                       temps_mean="mean temperature (lag 0)", 
                       ppts_warm= "warm precipitation (lag 0)", 
                       ppts_cool= "cool precipitation (lag 0)"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","k-rat accessible", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","k-rat accessible", "k-rat inaccessible"))+
  #ggtitle("male")+#geom_label(label="D", x=2, y=8, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con2=ggarrange(plot_con21, plot_con22+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                    font.label = list(size=12))
annotate_figure(plot_con2, top=text_grob("desert pocket mouse", face="bold"))

#biotic factors####

#PB plots####
plot_con31=dwplot(list(dipof_con3,pbf_con3, pbf_ex3),style="distribution",vars_order = c( "pbs","pps", "dipos"))%>%
  relabel_predictors(c( pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  ggtitle("female")+geom_label(label="A", x=5, y=6, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con32=dwplot(list(dipof_con3,pbm_con3, pbm_ex3),style="distribution",vars_order = c( "pbs","pps", "dipos"))%>%
  relabel_predictors(c( pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  ggtitle("male")+geom_label(label="B", x=5, y=5.7, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con3=ggarrange(plot_con31, plot_con32+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                    font.label = list(size=12))
annotate_figure(plot_con3, top=text_grob("Bailey's pocket mouse", face="bold"))

#PP plots####
plot_con41=dwplot(list(dipof_con3,ppf_con3, ppf_ex3),style="distribution",vars_order = c( "pbs","pps", "dipos"))%>%
  relabel_predictors(c( pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  ggtitle("female")+geom_label(label="C", x=4, y=5.5, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con42=dwplot(list(dipof_con3,ppm_con3, ppm_ex3),style="distribution",vars_order = c( "pbs","pps", "dipos"))%>%
  relabel_predictors(c( pbs="Bailey's pocket mouse biomass",pps= "desert pocket mouse biomass", dipos= "Dipodomys sp. biomass"))+
  scale_color_manual(name = "Species",
                     values=c("sky blue", "green", "orange"),
                     labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  scale_fill_manual(name = "Species",
                    values=c("sky blue", "green","orange"),
                    labels = c("kangaroo rats","control", "k-rat inaccessible"))+
  ggtitle("male")+geom_label(label="D", x=3, y=6, label.size=NA)+
  theme_classic()+geom_vline(xintercept=0, linetype="dotted")+theme_pubr(base_size = 10)

plot_con4=ggarrange(plot_con41, plot_con42+theme(axis.text.y = element_blank()),widths = c(0.35, 0.20), nrow=1, common.legend = T, legend="bottom",
                    font.label = list(size=12))
annotate_figure(plot_con4, top=text_grob("desert pocket mouse", face="bold"))
