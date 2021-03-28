####load packages####

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)

####load cleaned individual-level data####
Portal_data_indiv=summarize_individual_rodents(
  path = get_default_data_path(),
  clean = TRUE,
  type = "Rodents",
  length = "all",
  unknowns = FALSE,
  time = "date",
  fillweight = FALSE,
  min_plots = 1,
  min_traps = 1,
  download_if_missing = TRUE,
  quiet = FALSE
)%>%filter(!(treatment=="removal"), !is.na(treatment), !is.na(sex))

#PB DATASET####

#MALES####

portal_male=Portal_data_indiv%>%filter(sex=="M") 
head(portal_male)

repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

PB=repro_male%>%
  filter(species=="PB", wgt >=16)

#get count of reproductive males per month per year per trt
pb_dat=PB%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PB=portal_male%>%
  filter(species=="PB",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion=right_join(pb_dat, total_PB)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion[is.na(total_proportion)] <- 0 #set non-detects to 0

ggplot(total_proportion, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB males")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#FEMALE####

portal_female=Portal_data_indiv%>%filter(sex=="F") #49% of individuals are males
head(portal_male)

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | nipples==c("R", "E", "B") | lactation=="L")

PBf=repro_female%>%
  filter(species=="PB", wgt >=21)%>%
  arrange(wgt)

#get count of reproductive males per month per year per trt
pbf_dat=PBf%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_PBf=portal_female%>%
  filter(species=="PB",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_f=right_join(pbf_dat, total_PBf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_f[is.na(total_proportion_f)] <- 0 #set non-detects to 0

ggplot(total_proportion_f, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB females")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

PB_all=rbind(total_proportion, total_proportion_f)
PB_all=as.data.frame(PB_all)

#add ndvi monthly data####
prod=ndvi(level="monthly")

prod2=prod%>%
  mutate(year=year(date), month=month(date))

pb_all_prod=right_join(prod2,PB_all)

write.csv(pb_all_prod, "PB_MF_full_prod.csv")

#visualization####
ggplot(pb_all_prod, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("Bailey's pocket mouse")+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(PB_all, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB")+
#  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
