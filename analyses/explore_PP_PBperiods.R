#exploratory work on effect of temporal scale (presence of PB) on PP reprod######

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(ggpubr)

source("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/RScripts/data_cleaning_functions_Supp.R")

####load cleaned data####

Portal_reprod_indiv=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/Portal_reproductive_indiv.csv")

#pre-PB data####

Portal_prePB=Portal_reprod_indiv%>%
  filter(!(year>1995))

Portal_PB=Portal_reprod_indiv%>%
  filter(!(year<1996) , !(year>2010))

Portal_postPB=Portal_reprod_indiv%>%
  filter(!(year<2011) , !(year>2014))

#PP data####

#males####

#pre-PB period####

Portal_male_pre=Portal_prePB%>%filter(sex=="M", !is.na(sex), !is.na(treatment)) 
head(portal_male)

PP_repro_male_pre=Portal_male_pre%>%
  filter(testes==c("S", "M", "R"),species=="PP", wgt>=15)

#get count of reproductive males per month per year per trt
pp_dat_pre=PP_repro_male_pre%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PP_pre=Portal_male_pre%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_pre=right_join(pp_dat_pre, total_PP_pre)%>%
  mutate(proportion=reproductive/abundance, sex="male", period="pre")%>%
  arrange(proportion)

PP_proportion_pre[is.na(PP_proportion_pre)] <- 0 #set non-detects to 0

PPM_pre=ggplot(PP_proportion_pre, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#PB colonization period####

Portal_male_PB=Portal_PB%>%filter(sex=="M", !is.na(sex), !is.na(treatment)) 

PP_repro_male_PB=Portal_male_PB%>%
  filter(testes==c("S", "M", "R"),species=="PP", wgt>=13)

#get count of reproductive males per month per year per trt
pp_dat_PB=PP_repro_male_PB%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PP_PB=Portal_male_PB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_PB=right_join(pp_dat_PB, total_PP_PB)%>%
  mutate(proportion=reproductive/abundance, sex="male", period="during")%>%
  arrange(proportion)

PP_proportion_PB[is.na(PP_proportion_PB)] <- 0 #set non-detects to 0

PPm_PB=ggplot(PP_proportion_PB, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#post-PB period####
Portal_male_postPB=Portal_postPB%>%filter(sex=="M", !is.na(sex), !is.na(treatment)) 

PP_repro_male_postPB=Portal_male_postPB%>%
  filter(testes==c("S", "M", "R"),species=="PP", wgt>=14)

#get count of reproductive males per month per year per trt
pp_dat_postPB=PP_repro_male_postPB%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PP_postPB=Portal_male_postPB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_postPB=right_join(pp_dat_postPB, total_PP_postPB)%>%
  mutate(proportion=reproductive/abundance, sex="male", period="post")%>%
  arrange(proportion)

PP_proportion_postPB[is.na(PP_proportion_postPB)] <- 0 #set non-detects to 0

PPm_postPB=ggplot(PP_proportion_postPB, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#females####

#pre-PB period####

Portal_female_pre=Portal_prePB%>%filter(sex=="F", !is.na(sex), !is.na(treatment)) 

PP_repro_female_pre=Portal_female_pre%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L", species=="PP", wgt>=13)

#get count of reproductive males per month per year per trt
ppf_dat_pre=PP_repro_female_pre%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PPf_pre=Portal_female_pre%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_pre=right_join(ppf_dat_pre, total_PPf_pre)%>%
  mutate(proportion=reproductive/abundance, sex="female", period="pre")%>%
  arrange(proportion)

PPf_proportion_pre[is.na(PPf_proportion_pre)] <- 0 #set non-detects to 0

PPf_pre=ggplot(PPf_proportion_pre, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#PB colonization period####

Portal_female_PB=Portal_PB%>%filter(sex=="F", !is.na(sex), !is.na(treatment)) 

PP_repro_female_PB=Portal_female_PB%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L", species=="PP", wgt>=12)

#get count of reproductive males per month per year per trt
ppf_dat_PB=PP_repro_female_PB%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PPf_PB=Portal_female_PB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_PB=right_join(ppf_dat_PB, total_PPf_PB)%>%
  mutate(proportion=reproductive/abundance, sex="female", period="during")%>%
  arrange(proportion)

PPf_proportion_PB[is.na(PPf_proportion_PB)] <- 0 #set non-detects to 0

PPf_PB=ggplot(PPf_proportion_PB, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#post-PB period####
Portal_female_postPB=Portal_postPB%>%filter(sex=="F", !is.na(sex), !is.na(treatment)) 

PP_repro_female_postPB=Portal_female_postPB%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L", species=="PP", wgt>=12)

#get count of reproductive males per month per year per trt
ppf_dat_postPB=PP_repro_female_postPB%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PPf_postPB=Portal_female_postPB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_postPB=right_join(ppf_dat_postPB, total_PPf_postPB)%>%
  mutate(proportion=reproductive/abundance, sex="female", period="post")%>%
  arrange(proportion)

PPf_proportion_postPB[is.na(PPf_proportion_postPB)] <- 0 #set non-detects to 0

PPf_postPB=ggplot(PPf_proportion_postPB, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#combine all####
all_PP_periods=rbind(PPf_proportion_pre, PPf_proportion_PB, PPf_proportion_postPB,
                     PP_proportion_pre, PP_proportion_PB, PP_proportion_postPB)

#add covariates monthly data####

prod=ndvi(level="monthly", sensor="landsat", fill=TRUE)

prod2=prod%>%
  mutate(year=year(date), month=month(date))
#add lag to ndvi####
ndvi_lag=prod2%>%
  mutate(lag_ndvi=lag(ndvi,order_by=date))%>%
  filter(!(year<1988))

temp=weather(level="monthly", fill=TRUE, horizon=90)%>%
  select(year,month,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)

temp$date=as.Date(paste(temp$year, temp$month, 01), "%Y %m %d")

#add lag to temp####
temp_lag=temp%>%
  mutate(lag_temp_mean=lag(meantemp,order_by=date),
         lag_temp_min=lag(mintemp,order_by=date),
         lag_temp_max=lag(maxtemp,order_by=date),
         lag_ppt=lag(precipitation,order_by=date),
         lag_ppt_warm_2=lag(warm_precip,n=2,order_by=date),
         lag_ppt_cool_2=lag(cool_precip,n=2,order_by=date),
         lag_ppt_warm=lag(warm_precip,order_by=date),
         lag_ppt_cool=lag(cool_precip,order_by=date))%>%
  filter(!(year<1988))
#combine weather and PP data####

all_PP=right_join(ndvi_lag,all_PP_periods)
all_PP_temp=right_join(temp_lag, all_PP)

#add biomass data####
bmass=biomass(level="Plot", type="Rodents",
              clean=TRUE, plots="all", time="date", shape="crosstab")%>%
  mutate(month=month(censusdate), date=day(censusdate), year=year(censusdate))%>%
  filter(!(year<1988), !(year>2014), !(treatment %in%c("spectabs", "removal")))

DO_bmass=bmass%>%select(DO, plot, treatment, month, year, date)%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_DO=sum(DO))

DS_bmass=bmass%>%select(DS, plot, treatment, month, year, date)%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_DS=sum(DS))

DM_bmass=bmass%>%select(DM, plot, treatment, month, year, date)%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_DM=sum(DM))

PP_bmass=bmass%>%select(PP, plot, treatment, month, year, date)%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_PP=sum(PP))

PB_bmass=bmass%>%select(PB, plot, treatment, month, year, date)%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_PB=sum(PB))

all_bmass11=left_join(all_PP_temp,DM_bmass, by=c("month", "year", "treatment"))
all_bmass12=left_join(all_bmass11,DO_bmass, by=c("month", "year", "treatment"))
all_bmass13=left_join(all_bmass12,DS_bmass, by=c("month", "year", "treatment"))
all_bmass2=left_join(all_bmass13,PB_bmass, by=c("month", "year", "treatment"))
all_bmass3=left_join(all_bmass2,PP_bmass, by=c("month", "year", "treatment"))%>%
  mutate(DIPO_bmass=rowSums(.[26:28]))%>%drop_na()

all_bmass3$years=(all_bmass3$year-mean(all_bmass3$year))/(2*sd(all_bmass3$year))
all_bmass3$ndvis=(all_bmass3$ndvi-mean(all_bmass3$ndvi))/(sd(all_bmass3$ndvi))
all_bmass3$ndvis_lag=(all_bmass3$lag_ndvi-mean(all_bmass3$lag_ndvi))/(sd(all_bmass3$lag_ndvi))
all_bmass3$temps_lag_mean=(all_bmass3$lag_temp_mean-mean(all_bmass3$lag_temp_mean))/(2*sd(all_bmass3$lag_temp_mean))
all_bmass3$temps_lag_min=(all_bmass3$lag_temp_min-mean(all_bmass3$lag_temp_min))/(2*sd(all_bmass3$lag_temp_min))
all_bmass3$temps_lag_max=(all_bmass3$lag_temp_max-mean(all_bmass3$lag_temp_max))/(2*sd(all_bmass3$lag_temp_max))
all_bmass3$temps_mean=(all_bmass3$meantemp-mean(all_bmass3$meantemp))/(2*sd(all_bmass3$meantemp))
all_bmass3$temps_min=(all_bmass3$mintemp-mean(all_bmass3$mintemp))/(2*sd(all_bmass3$mintemp))
all_bmass3$temps_max=(all_bmass3$maxtemp-mean(all_bmass3$maxtemp))/(2*sd(all_bmass3$maxtemp))
all_bmass3$ppts=(all_bmass3$precipitation-mean(all_bmass3$precipitation))/(2*sd(all_bmass3$precipitation))
all_bmass3$ppts_lag=(all_bmass3$lag_ppt-mean(all_bmass3$lag_ppt))/(2*sd(all_bmass3$lag_ppt))
all_bmass3$dipos=(all_bmass3$DIPO_bmass-mean(all_bmass3$DIPO_bmass))/(2*sd(all_bmass3$DIPO_bmass))
all_bmass3$pbs=(all_bmass3$bmass_PB-mean(all_bmass3$bmass_PB))/(2*sd(all_bmass3$bmass_PB))
all_bmass3$pps=(all_bmass3$bmass_PP-mean(all_bmass3$bmass_PP))/(2*sd(all_bmass3$bmass_PP))
all_bmass3$ppts_warm=(all_bmass3$warm_precip-mean(all_bmass3$warm_precip))/(2*sd(all_bmass3$warm_precip))
all_bmass3$ppts_cool=(all_bmass3$cool_precip-mean(all_bmass3$cool_precip))/(2*sd(all_bmass3$cool_precip))
all_bmass3$ppts_lag_warm=(all_bmass3$lag_ppt_warm-mean(all_bmass3$lag_ppt_warm))/(2*sd(all_bmass3$lag_ppt_warm))
all_bmass3$ppts_lag_cool=(all_bmass3$lag_ppt_cool-mean(all_bmass3$lag_ppt_cool))/(2*sd(all_bmass3$lag_ppt_cool))
all_bmass3$ppts_lag_warm_2=(all_bmass3$lag_ppt_warm_2-mean(all_bmass3$lag_ppt_warm_2))/(2*sd(all_bmass3$lag_ppt_warm_2))
all_bmass3$ppts_lag_cool_2=(all_bmass3$lag_ppt_cool_2-mean(all_bmass3$lag_ppt_cool_2))/(2*sd(all_bmass3$lag_ppt_cool_2))

#period and sex-specific####
ppf_con_pre=all_bmass3%>%filter(sex=="female", period=="pre", treatment=="control")
ppm_con_pre=all_bmass3%>%filter(sex=="male", period=="pre", treatment=="control")


#phenology visualization####
PP_pre=ggarrange(PPM_pre,PPf_pre, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))
annotate_figure(PP_pre, fig.lab="pre-PB period (1977-1994)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")


PP_PB=ggarrange(PPm_PB,PPf_PB, nrow=1, common.legend = T, legend="right",
                 font.label = list(size=12))
annotate_figure(PP_PB, fig.lab="PB colonization period (1995-2010)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")

PP_postPB=ggarrange(PPm_postPB,PPf_postPB, nrow=1, common.legend = T, legend="right",
                font.label = list(size=12))
annotate_figure(PP_postPB, fig.lab="post-PB period 2011-2014)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")

#MODELS####

#control####
ppf_con_pre=all_bmass3%>%filter(sex=="female", period=="pre", treatment=="control")
ppm_con_pre=all_bmass3%>%filter(sex=="male", period=="pre", treatment=="control")

ppf_con_mid=all_bmass3%>%filter(sex=="female", period=="during", treatment=="control")
ppm_con_mid=all_bmass3%>%filter(sex=="male", period=="during", treatment=="control")

ppf_con_post=all_bmass3%>%filter(sex=="female", period=="post", treatment=="control")
ppm_con_post=all_bmass3%>%filter(sex=="male", period=="post", treatment=="control")

#exclosure####
ppf_ex_pre=all_bmass3%>%filter(sex=="female", period=="pre", treatment=="exclosure")
ppm_ex_pre=all_bmass3%>%filter(sex=="male", period=="pre", treatment=="exclosure")

ppf_ex_mid=all_bmass3%>%filter(sex=="female", period=="during", treatment=="exclosure")
ppm_ex_mid=all_bmass3%>%filter(sex=="male", period=="during", treatment=="exclosure")

ppf_ex_post=all_bmass3%>%filter(sex=="female", period=="post", treatment=="exclosure")
ppm_ex_post=all_bmass3%>%filter(sex=="male", period=="post", treatment=="exclosure")

ppm_ex1_mid=mgcv::gam(proportion~s(month,bs="cc")+year+ndvis+ndvis_lag+temps_mean+
                        temps_lag_mean+ppts_warm+ppts_lag_warm+ppts_cool+ppts_lag_cool, 
                      data=ppm_ex_mid, method = 'REML', 
                      weights = abundance, family = binomial)
summary(ppm_ex1_mid)

plot(ppf_ex_pre$proportion~ppf_ex_pre$year)

