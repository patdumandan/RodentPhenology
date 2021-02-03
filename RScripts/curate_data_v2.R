library(dplyr)
library(tidyr)
library(portalr)

Portal_data=load_rodent_data()
Portal_rodent=Portal_data[["rodent_data"]]
Portal_plots=Portal_data[["plots_table"]]

portal1=inner_join(Portal_rodent, Portal_plots, by= c("plot", "month", "year"))%>%
  filter(!(treatment=="removal"))%>%
  mutate(month=factor(month), Month=recode(month, "1"= "Jan", "2"="Feb", "3"="Mar", "4"="Apr",
                                           "5"="May","6"="Jun", "7"="Jul", "8"="Aug", "9"="Sept",
                                           "10"="Oct","11"="Nov", "12"="Dec"))%>%
  select(period, month, Month, day, year, plot, treatment, species, sex, reprod, vagina, nipples,lactation, pregnant, testes,hfl,wgt)
str(portal1)
unique(portal1$species)
unique(portal1$treatment)

#male individuals####
#spetabs an removal plots--what to do??

target_repro=c("S", "M", "R")

#size thresholds####
BA=portal1%>%
  filter(species=="BA",testes %in% target_repro, wgt >=6)%>%
  arrange(wgt)

DM=portal1%>%
  filter(species=="DM",testes %in% target_repro, wgt >=15)%>%
  arrange(wgt)

DO=portal1%>%
  filter(species=="DO",testes %in% target_repro, wgt >=29)%>%
  arrange(wgt)

DS=portal1%>%
  filter(species=="DS",testes %in% target_repro, wgt >=12)%>%
  arrange(wgt)

NEA=portal1%>%
  filter(species=="NA",testes %in% target_repro, wgt >=121)%>%
  arrange(wgt)

OL=portal1%>%
  filter(species=="OL",testes %in% target_repro, wgt >=19)%>%
  arrange(wgt)

OT=portal1%>%
  filter(species=="OT",testes %in% target_repro, wgt >=10)%>%
  arrange(wgt)

PH=portal1%>%
  filter(species=="PH",testes %in% target_repro, wgt >=18)%>%
  arrange(wgt)

PL=portal1%>%
  filter(species=="PL",testes %in% target_repro, wgt >=20)%>%
  arrange(wgt)

PB=portal1%>%
  filter(species=="PB",testes %in% target_repro,wgt >=16)%>%
  arrange(wgt)

PP=portal1%>%
  filter(species=="PP",testes %in% target_repro, wgt >=10)%>%
  arrange(wgt)

PE=portal1%>%
  filter(species=="PE",testes %in% target_repro, wgt >=7)%>%
  arrange(wgt)

PF=portal1%>%
  filter(species=="PF",testes %in% target_repro, wgt >=4)%>%
  arrange(wgt)

PM=portal1%>%
  filter(species=="PM",testes %in% target_repro, wgt >=11)%>%
  arrange(wgt)
#1 individual identified as juv but size within reproductive??

RM=portal1%>%
  filter(species=="RM",testes %in% target_repro, wgt >=4)%>%
  arrange(wgt)

RF=portal1%>%
  filter(species=="RF",testes %in% target_repro, wgt >=11)%>%
  arrange(wgt)

RO=portal1%>%
  filter(species=="RO",testes %in% target_repro, wgt >=6)%>%
  arrange(wgt)

SF=portal1%>%
  filter(species=="SF",testes %in% target_repro, wgt >=39)%>%
  arrange(wgt)

SH=portal1%>%
  filter(species=="SH",testes %in% target_repro, wgt >=51)%>%
  arrange(wgt)

SO=portal1%>%
  filter(species=="SO",testes %in% target_repro, wgt >=68)%>%
  arrange(wgt)

#create full dataset with all species

full_repro_male_dat=rbind(SO,SH,SF,RO,RM,RF,PP,PM,PL,PH,PF,PE,PB,OT,OL,NEA,DS,DO,DM,BA)
write.csv(full_repro_male_dat, "repro_male_all.csv")

repro_male_all=read.csv("repro_male_all.csv")

#read in individual rodents
indivs=summarize_individual_rodents()%>%
  filter(!(treatment=="removal"& NA), sex=="M")

#summarize repro_male count
repro_dat=repro_male_all%>%
  group_by(month, year, treatment, species)%>%
  summarise(reproductive=n())

#summarize total rodents count
total_rodents=indivs%>%
  group_by(month, year, treatment, species)%>%
  summarise(abundance=n())

total_proportion=left_join(repro_dat, total_rodents)%>%
  mutate(proportion=reproductive/abundance)

#figure out how to address the NA in species

#bad periods???

