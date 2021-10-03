#exploratory work on effect of temporal scale (presence of PB) on PP reprod######

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)

source("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/RScripts/data_cleaning_functions_Supp.R")

####load cleaned data####

Portal_reprod_indiv=read.csv("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/ReproData/Portal_reproductive_indiv.csv")

#pre-PB data####

Portal_prePB=Portal_reprod_indiv%>%
  filter(!(year>1997))

Portal_PB=Portal_reprod_indiv%>%
  filter(!(year<1998) , !(year>2010))

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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PP_pre=Portal_male_pre%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_pre=right_join(pp_dat_pre, total_PP_pre)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PP_PB=Portal_male_PB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_PB=right_join(pp_dat_PB, total_PP_PB)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PP_postPB=Portal_male_postPB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PP_proportion_postPB=right_join(pp_dat_postPB, total_PP_postPB)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PPf_pre=Portal_female_pre%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_pre=right_join(ppf_dat_pre, total_PPf_pre)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PPf_PB=Portal_female_PB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_PB=right_join(ppf_dat_PB, total_PPf_PB)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
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
  summarise(reproductive=n(month))

#get total observed abundance for MALES per month per year per trt
total_PPf_postPB=Portal_female_postPB%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n(month))

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
PPf_proportion_postPB=right_join(ppf_dat_postPB, total_PPf_postPB)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
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

#phenology visualization####
PP_pre=ggarrange(PPM_pre,PPf_pre, nrow=1, common.legend = T, legend="right",
              font.label = list(size=12))
annotate_figure(PP_pre, fig.lab="pre-PB period (1977-1997)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")


PP_PB=ggarrange(PPm_PB,PPf_PB, nrow=1, common.legend = T, legend="right",
                 font.label = list(size=12))
annotate_figure(PP_PB, fig.lab="PB colonization period (1998-2010)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")

PP_postPB=ggarrange(PPm_postPB,PPf_postPB, nrow=1, common.legend = T, legend="right",
                font.label = list(size=12))
annotate_figure(PP_postPB, fig.lab="post-PB period 2011-2014)", fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab.face = "bold")


