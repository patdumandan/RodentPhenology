####load packages####

library(dplyr)
library(tidyr)
library(portalr)
library(ggplot2)
library(lubridate)
library(reshape2)

source("https://raw.githubusercontent.com/patdumandan/ReproPhenology/main/RScripts/data_cleaning_functions_Supp.R")
####load cleaned data####
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

Portal_rodent=read.csv("https://raw.githubusercontent.com/weecology/PortalData/main/Rodents/Portal_rodent.csv")
Portal_data_indiv=left_join(Portal_data, Portal_rodent)%>%
  select(period, month, day, year, treatment, plot, stake, species, sex, reprod, age, testes, vagina, pregnant, nipples, lactation,
         hfl, wgt,tag,note2, note5)

#assign tag IDs for untagged individuals (0 and NA in tag column)####

all_tag=id_unknowns(Portal_data_indiv, 19) 

#find and remove bad periods (periods with only one day of trapping)####

Portal_rodent_trapping= read.csv("https://raw.githubusercontent.com/weecology/PortalData/main/Rodents/Portal_rodent_trapping.csv")
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

#write.csv(Portal_clean, "Portal_reproductive_indiv.csv")

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
  filter(species=="DM", wgt>=27)

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
  mutate(species="DM")%>%
  filter(treatment=="control")

#DO dataset####
#males####
repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

DO=repro_male%>%
  filter(species=="DO", wgt>=34)

#get count of reproductive males per month per year per trt
do_dat=DO%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_DO=portal_male%>%
  filter(species=="DO")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_do_m=right_join(do_dat, total_DO)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_do_m[is.na(total_proportion_do_m)] <- 0 #set non-detects to 0

#females####

portal_female=Portal_clean%>%filter(sex=="F") #49% of individuals are males

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L")

DOf=repro_female%>%
  filter(species=="DO", wgt >=28)

#get count of reproductive males per month per year per trt
dof_dat=DOf%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_DOf=portal_female%>%
  filter(species=="DO")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_do_f=right_join(dof_dat, total_DOf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_do_f[is.na(total_proportion_do_f)] <- 0 #set non-detects to 0

DO_all=rbind(total_proportion_do_m, total_proportion_f)
DO_all=as.data.frame(DO_all)%>%
  mutate(species="DO")%>%filter(treatment=="control")

#DS dataset####
#males####
repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))

DS=repro_male%>%
  filter(species=="DS", wgt>=100)

#get count of reproductive males per month per year per trt
ds_dat=DS%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for MALES per month per year per trt
total_DS=portal_male%>%
  filter(species=="DS")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_ds_m=right_join(ds_dat, total_DS)%>%
  mutate(proportion=reproductive/abundance, sex="male")%>%
  arrange(proportion)

total_proportion_ds_m[is.na(total_proportion_ds_m)] <- 0 #set non-detects to 0

#females####

portal_female=Portal_clean%>%filter(sex=="F") #49% of individuals are males

repro_female=portal_female%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | 
           nipples==c("R", "E", "B") | lactation=="L")

DSf=repro_female%>%
  filter(species=="DS", wgt >=78)

#get count of reproductive males per month per year per trt
dsf_dat=DSf%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_DSf=portal_female%>%
  filter(species=="DS")%>%
  group_by(month,year, treatment)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion_ds_f=right_join(dsf_dat, total_DSf)%>%
  mutate(proportion=reproductive/abundance, sex="female")%>%
  arrange(proportion)

total_proportion_ds_f[is.na(total_proportion_ds_f)] <- 0 #set non-detects to 0

DS_all=rbind(total_proportion_ds_m, total_proportion_ds_f)
DS_all=as.data.frame(DS_all)%>%
  mutate(species="DS")%>%filter(treatment=="control")

#combine PB, PP and DM datasets####
all_sp=rbind(PB_all, PP_all, DM_all, DO_all, DS_all)

#write.csv(all_sp, "PB_PP_dipos_reprod.csv")

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

#combine weather and rodent data####
all_prod=right_join(ndvi_lag,all_sp)
all_prod_temp=right_join(temp_lag, all_prod)

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

all_bmass11=left_join(all_prod_temp,DM_bmass, by=c("month", "year", "treatment"))
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

#write.csv(all_bmass3, "reproductive_full_data.csv")

#visualization####
pb_plot=all_bmass3%>%filter(species=="PB")
pp_plot=all_bmass3%>%filter(species=="PP")
dm_plot=all_bmass3%>%filter(species=="DM")
do_plot=all_bmass3%>%filter(species=="DO")
ds_plot=all_bmass3%>%filter(species=="DS")

#PP and PB####
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

#Dipos####
ggplot(DM_all, aes(y=proportion, x=month, col=sex)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DM")+ 
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(DO_all, aes(y=proportion, x=month, col=sex)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DO")+ 
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(DS_all, aes(y=proportion, x=month, col=sex)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("DS")+ 
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#template code to obtain data for species and sex-specific models####
PP_male_con=pp_plot%>%filter(treatment=="control", sex=="male")
PP_male_ex=pp_plot%>%filter(treatment=="exclosure", sex=="male")
PP_female_con=pp_plot%>%filter(treatment=="control", sex=="female")
PP_female_ex=pp_plot%>%filter(treatment=="exclosure", sex=="female")