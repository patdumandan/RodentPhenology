#try it with female PBs
require(dplyr)
tot_dat=read.csv("./reconfigured_data/raw_cleaned/portal_reproductive.csv")
fem_dat=tot_dat%>%filter(sex=="F") #51% of all repros are female


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
head(total_proportion)
#length(unique(total_proportion$species)) #21 spp
#max(total_proportion$proportion, na.rm=T) #1