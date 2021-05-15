####load packages####

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)

source("./RScripts/data_cleaning_functions_Supp.R")

####load cleaned individual-level data####
Portal_data=summarize_individual_rodents(
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
#Note: 18 plots included based on Ellen's paper

#add note5 column to filter out dead indivs####

Portal_rodent=read.csv("./PortalData/Rodents/Portal_rodent.csv")
Portal_data_indiv=left_join(Portal_data, Portal_rodent)%>%
  select(period, month, day, year, treatment, plot, stake, species, sex, reprod, age, testes, vagina, pregnant, nipples, lactation,
         hfl, wgt,tag,note2, note5)

#assign tag IDs for untagged individuals (0 and NA in tag column)####

all_tag=id_unknowns(Portal_data_indiv, 19) 

#find and remove bad periods (periods with only one day of trapping)####

Portal_rodent_trapping= read.csv("./PortalData/Rodents/Portal_rodent_trapping.csv")
tdat=Portal_rodent_trapping%>%group_by(period)%>%summarise(count=sum(sampled))%>%arrange(count)
bad_periods <- filter(tdat, count < 20) #based on Sarah's code
bad_periods <- as.list(bad_periods$period)

Portal_no_badperiod=all_tag[-which(all_tag$period %in%bad_periods),]%>%
  mutate(tag=as.character(tag)) #necessary so function for starred and duplicate tags will not break

#check quality of tags####

#make sure that records with * in note2 are recognized as new individuals

tags = unique(Portal_no_badperiod$tag)
star_tags = starred_tags(dat=Portal_no_badperiod, tags=tags, spp_col=8, tag_col=19)

#locate dead individuals
tags=unique(star_tags$tag)
dead_dat=is_dead(dat=star_tags, tags=tags, spp_col=8, tag_col=19)

# locate records of duplicated tags
tags=unique(dead_dat$tag)
dup_dat= is_duplicate_tag(dat=dead_dat, tags=tags, spp_col=8, tag_col=19) #returns a list of 2
duptags = unique(dup_dat$bad$tag)
no_dup = dup_dat$data[-which(dup_dat$data$tag %in% duptags),] #delete rows flagged as duplicates without clear resolution

# identify records where multiple indivs share same tag in same period
tags = unique(no_dup$tag)
same = same_period(no_dup, tags)

#eliminate tags that appear more than once in the same period - questionable data
sametags = unique(same$tag)
Portal_no_same= no_dup[-which(no_dup$tag %in% sametags),]

# "clean" data

Portal_clean=subsetDat(Portal_no_same)

#Note: this analysis does not necessarily follow the capture history of 
#individuals, what we want are just the event IDs/observations of reprod.
#characteristics to determine peak timing of breeding events

#PB DATASET####

#males####

portal_male=Portal_clean%>%filter(sex=="M", !is.na(sex), !is.na(treatment)) 
head(portal_male)

repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

PB=repro_male%>%
  filter(species=="PB", wgt>=18)

#get count of reproductive males per month per year per trt
pb_dat=PB%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PB=portal_male%>%
  filter(species=="PB")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion=right_join(pb_dat, total_PB)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion[is.na(total_proportion)] <- 0 #set non-detects to 0

#females####

portal_female=Portal_clean%>%filter(sex=="F") #49% of individuals are males
head(portal_male)

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L")

PBf=repro_female%>%
  filter(species=="PB", wgt >=21)%>%
  arrange(wgt)

#get count of reproductive males per month per year per trt
pbf_dat=PBf%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_PBf=portal_female%>%
  filter(species=="PB")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_f=right_join(pbf_dat, total_PBf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_f[is.na(total_proportion_f)] <- 0 #set non-detects to 0

PB_all=rbind(total_proportion, total_proportion_f)
PB_all=as.data.frame(PB_all)%>%
  mutate(species="PB")

#PP DATASET####

#males####

PP=repro_male%>%
  filter(species=="PP", wgt>=13)%>%
  arrange(wgt)

#get count of reproductive males per month per year per trt
PP_dat=PP%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_PP=portal_male%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male were recorded
total_proportion_pp_m=right_join(PP_dat, total_PP)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_pp_m[is.na(total_proportion_pp_m)] <- 0 #set non-detects to 0

#females####

portal_female=Portal_clean%>%filter(sex=="F") #49% of individuals are males
head(portal_male)

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | nipples==c("R", "E", "B") | lactation=="L")

PPf=repro_female%>%
  filter(species=="PP",  wgt >=12)%>%
  arrange(wgt)

#get count of reproductive males per month per year per trt
PPf_dat=PPf%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_PPf=portal_female%>%
  filter(species=="PP")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_pp_f=right_join(PPf_dat, total_PPf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_pp_f[is.na(total_proportion_pp_f)] <- 0 #set non-detects to 0

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
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_DM=portal_male%>%
  filter(species=="DM")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_dm_m=right_join(DM_dat, total_DM)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_dm_m[is.na(total_proportion_dm_m)] <- 0 #set non-detects to 0

#females####

DMf=repro_female%>%
  filter(species=="DM", wgt>=27)%>%arrange(wgt)

#get count of reproductive males per month per year per trt
DMf_dat=DMf%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_DMf=portal_female%>%
  filter(species=="DM")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_dm_f=right_join(DMf_dat, total_DMf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_dm_f[is.na(total_proportion_dm_f)] <- 0 #set non-detects to 0

DM_all=rbind(total_proportion_dm_m, total_proportion_dm_f)
DM_all=as.data.frame(DM_all)%>%
  mutate(species="DM")

#combine PB, PP and DM datasets####
all_sp=rbind(PB_all, PP_all, DM_all)

#add ndvi and ppt monthly data####

prod=ndvi(level="monthly", sensor="landsat", fill=TRUE)

prod2=prod%>%
  mutate(year=year(date), month=month(date))

ppt=weather(level="monthly", fill=TRUE)%>%select(year,month,precipitation)

all_prod=right_join(prod2,all_sp)
all_prod_ppt=right_join(ppt, all_prod)%>%filter(!is.na(precipitation), !is.na(ndvi))

#add biomass data####
bmass=biomass(level="Plot", type="Rodents",
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

all_bmass1=left_join(all_prod_ppt,DM_bmass, by=c("month", "year", "treatment"))
all_bmass2=left_join(all_bmass1,PB_bmass, by=c("month", "year", "treatment"))
all_bmass3=left_join(all_bmass2,PP_bmass, by=c("month", "year", "treatment"))

# adding lags of weather variables####
var_lag=all_bmass3%>%
  mutate(month=as.integer(month), lag_month=month-1)%>%
  select(lag_month, ndvi, precipitation, month, year)%>%
  rename(lag_ndvi=ndvi, lag_ppt=precipitation)

var_lag2=left_join(all_bmass3, var_lag, by=c("month"="lag_month", "year"))%>%
  distinct()

#visualization####
pb_plot=all3_bmass%>%filter(species=="PB")
pp_plot=all3_bmass%>%filter(species=="PP")
dm_plot=all3_bmass%>%filter(species=="DM")


ggplot(pb_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB")+ facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(pp_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PP")+facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(dm_plot, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DM")+facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
