####load packages####

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)

####load cleaned individual-level data####
Portal_data_indiv=summarize_individual_rodents(
  clean = TRUE,
  type = "Rodents",
  length = "all",
  unknowns = FALSE,
  fillweight = FALSE,
  min_plots = 1,
  min_traps = 1,
  download_if_missing = TRUE,
  quiet = FALSE
)%>%filter(!is.na(sex),!(treatment=="spectabs"), !(year<1988), !(year>2014), 
           plot %in%c(1, 2, 4, 8, 9, 11, 12, 14, 17, 22,3, 6, 13, 15, 18, 19, 20, 21))

#find and remove bad periods (periods with only one day of trapping)###
tdat=Portal_rodent_trapping%>%group_by(period)%>%summarise(count=sum(sampled))%>%arrange(count)
bad_periods <- filter(tdat, count < 20)
bad_periods <- as.list(bad_periods$period)

Portal_no_badperiod=Portal_data_indiv[-which(Portal_data_indiv$period %in%bad_periods),]

#PB DATASET####

#males####

portal_male=Portal_no_badperiod%>%filter(sex=="M") 
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

#females####

portal_female=Portal_no_badperiod%>%filter(sex=="F") #49% of individuals are males
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
PB_all=as.data.frame(PB_all)%>%
  mutate(species="PB")

#PP DATASET####

#males####

portal_male=Portal_no_badperiod%>%filter(sex=="M") 
head(portal_male)

repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

PP=repro_male%>%
  filter(species=="PP", wgt>=10)

#get count of reproductive males per month per year per trt
PP_dat=PP%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PP=portal_male%>%
  filter(species=="PP",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_pp_m=right_join(PP_dat, total_PP)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_pp_m[is.na(total_proportion_pp_m)] <- 0 #set non-detects to 0

ggplot(total_proportion_pp_m, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PP males")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#females####

portal_female=Portal_no_badperiod%>%filter(sex=="F") #49% of individuals are males
head(portal_male)

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | nipples==c("R", "E", "B") | lactation=="L")

PPf=repro_female%>%
  filter(species=="PP", wgt >=11)

#get count of reproductive males per month per year per trt
PPf_dat=PPf%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_PPf=portal_female%>%
  filter(species=="PP",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_pp_f=right_join(PPf_dat, total_PPf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_pp_f[is.na(total_proportion_pp_f)] <- 0 #set non-detects to 0

ggplot(total_proportion_pp_f, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PP females")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

PP_all=rbind(total_proportion_pp_m, total_proportion_pp_f)
PP_all=as.data.frame(PP_all)%>%
  mutate(species="PP")

#DM DATASET####

#males####

portal_male=Portal_no_badperiod%>%filter(sex=="M") 
head(portal_male)

repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

DM=repro_male%>%
  filter(species=="DM", wgt>=22)

#get count of reproductive males per month per year per trt
DM_dat=DM%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_DM=portal_male%>%
  filter(species=="DM",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_dm_m=right_join(DM_dat, total_DM)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_dm_m[is.na(total_proportion_dm_m)] <- 0 #set non-detects to 0

ggplot(total_proportion_dm_m, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DM males")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#females####

portal_female=Portal_no_badperiod%>%filter(sex=="F") #49% of individuals are males
head(portal_male)

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | nipples==c("R", "E", "B") | lactation=="L")

DMf=repro_female%>%
  filter(species=="DM", wgt>=25)

#get count of reproductive males per month per year per trt
DMf_dat=DMf%>%
  filter(!(treatment=="spectabs"))%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_DMf=portal_female%>%
  filter(species=="DM",!(treatment=="spectabs"))%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_dm_f=right_join(DMf_dat, total_DMf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_dm_f[is.na(total_proportion_dm_f)] <- 0 #set non-detects to 0

ggplot(total_proportion_dm_f, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DM females")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DM_all=rbind(total_proportion_dm_m, total_proportion_dm_f)
DM_all=as.data.frame(DM_all)%>%
  mutate(species="DM")

#combine PB, PP and DM datasets####
all_sp=rbind(PB_all, PP_all, DM_all)

#add ndvi and ppt monthly data####

prod=ndvi(level="monthly", sensor="gimms", fill=TRUE)

prod2=prod%>%
  mutate(year=year(date), month=month(date))
ppt=weather(level="monthly", fill=TRUE)%>%select(year,month,precipitation)

all_prod=right_join(prod2,all_sp)
all_prod_ppt=right_join(ppt, all_prod)%>%filter(!is.na(precipitation), !is.na(ndvi))
#write.csv(all_prod_ppt, "all_full_prod.csv")

#add biomass data####
bmass=biomass(path="repo", level="Plot", type="Rodents",
              clean=TRUE, plots="all", time="date", shape="crosstab")


DM_bmass=bmass%>%select(DM, plot, treatment, censusdate)%>%
  filter(!(treatment%in%c("removal", "spectabs")))%>%
  mutate(month=month(censusdate), date=day(censusdate), year=year(censusdate))%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_DM=sum(DM))

PB_bmass=bmass%>%select(PB, plot, treatment, censusdate)%>%
  filter(!(treatment%in%c("removal", "spectabs")))%>%
  mutate(month=month(censusdate), date=day(censusdate), year=year(censusdate))%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_PB=sum(PB))

PP_bmass=bmass%>%select(PP, plot, treatment, censusdate)%>%
  filter(!(treatment%in%c("removal", "spectabs")))%>%
  mutate(month=month(censusdate), date=day(censusdate), year=year(censusdate))%>%
  group_by(month, year, treatment)%>%
  summarise(bmass_PP=sum(PP))

all_bmass=left_join(all_prod_ppt, DM_bmass, by=c("month", "year", "treatment"))
all2_bmass=left_join(all_bmass, PB_bmass, by=c("month", "year", "treatment"))
all3_bmass=left_join(all2_bmass, PP_bmass, by=c("month", "year", "treatment"))%>%
  filter(!is.na(bmass_DM), !is.na(bmass_PB), !is.na(bmass_PP))


#visualization####
pb_plot=all3_bmass%>%filter(species=="PB")
pp_plot=all3_bmass%>%filter(species=="PP")
dm_plot=all3_bmass%>%filter(species=="DM")

ggplot(pb_plot, aes(y=proportion, x=month, col=treatment)) +
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

ggplot(pp_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("desert pocket mouse")+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(dm_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("Merriam's kangaroo rat")+
  facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

