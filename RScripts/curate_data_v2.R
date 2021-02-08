####load packages####

library(dplyr)
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
)%>%filter(!(treatment=="removal")& !is.na(treatment)& !is.na(sex))%>%
  mutate(month=as.character(month))

# add columns for reproductive traits
Portal_data=load_rodent_data(clean = TRUE)
Portal_rodent=Portal_data[["rodent_data"]]%>%mutate(month=factor(month))%>%
  filter(!is.na(sex))

#create full dataset without removal plots
portal1=left_join(Portal_data_indiv, Portal_rodent)%>%
  filter(!(treatment=="removal")& !is.na(treatment)& !is.na(sex))%>%
  mutate(month=as.character(month), Month=recode(month, "1"= "Jan", "2"="Feb", "3"="Mar", "4"="Apr",
                                           "5"="May","6"="Jun", "7"="Jul", "8"="Aug", "9"="Sept",
                                           "10"="Oct","11"="Nov", "12"="Dec"))%>%
  select(period, month, Month, day, year, plot, stake,
         treatment, species, sex, reprod, vagina, nipples,lactation, pregnant, testes,hfl,wgt, tag)

str(portal1)
#length(unique(portal1$species))
#unique(portal1$treatment)


#CREATE DATASET FOR MALE ADULTS ONLY

portal_male=portal1%>%filter(sex=="M") #49% of individuals are males
head(portal_male)

####determine threshold for breeding adult male individuals####

target_repro=c("S", "M", "R")
repro_male=portal_male%>%
  filter(testes==c("S", "M", "R"))
head(repro_male)

#BA=repro_male%>%
#  filter(species=="BA")%>%
#  arrange(wgt)

#size thresholds####
BA=repro_male%>%
  filter(species=="BA", wgt >=6)

DM=repro_male%>%
  filter(species=="DM", wgt >=15)

DO=repro_male%>%
  filter(species=="DO", wgt >=29)

DS=repro_male%>%
  filter(species=="DS", wgt >=12)

NEA=repro_male%>%
  filter(species=="NA", wgt >=121)

OL=repro_male%>%
  filter(species=="OL", wgt >=19)

OT=repro_male%>%
  filter(species=="OT", wgt >=10)

PH=repro_male%>%
  filter(species=="PH", wgt >=18)

PL=repro_male%>%
  filter(species=="PL", wgt >=20)

PB=repro_male%>%
  filter(species=="PB", wgt >=16)

PP=repro_male%>%
  filter(species=="PP", wgt >=10)

PE=repro_male%>%
  filter(species=="PE", wgt >=7)

PI=repro_male%>%
  filter(species=="PI", wgt >=15)

PF=repro_male%>%
  filter(species=="PF", wgt >=4)

PM=repro_male%>%
  filter(species=="PM", wgt >=11)

RM=repro_male%>%
  filter(species=="RM", wgt >=4)

RF=repro_male%>%
  filter(species=="RF", wgt >=11)

RO=repro_male%>%
  filter(species=="RO", wgt >=6)

SF=repro_male%>%
  filter(species=="SF", wgt >=39)

SH=repro_male%>%
  filter(species=="SH", wgt >=51)

SO=repro_male%>%
  filter(species=="SO", wgt >=68)

#create full dataset with all species

full_repro_male_dat=rbind(SO,SH,SF,RO,RM,RF,PP,PM,PL,PH,PF,PE,PB,OT,OL,NEA,DS,DO,DM,BA, PI)
full_repro_male_dat=as.data.frame(full_repro_male_dat)
head(full_repro_male_dat)

####calculate proportion of reproductive individuals####

#get count of reproductive males for each species per month per year per trt
repro_dat=full_repro_male_dat%>%
  group_by(month, year, treatment, species)%>%
  summarise(reproductive=n())

#get total observed abundance for each species per month per year per trt
total_rodents=portal_male%>%
  group_by(month,year, treatment, species)%>%
  summarise(abundance=n())

#calculate proportion
#this creates NAs for months when no reproductive male was recorded
total_proportion=right_join(repro_dat, total_rodents)%>%
  mutate(proportion=reproductive/abundance)%>%
  arrange(proportion)
length(unique(total_proportion$species))
max(total_proportion$proportion, na.rm=T)


