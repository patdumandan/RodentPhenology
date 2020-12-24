library(dplyr)
Portal_rodent=read.csv("./PortalData/Rodents/Portal_rodent.csv")
Portal_plots=read.csv("./PortalData/SiteandMethods/Portal_plots.csv")
portal1=inner_join(Portal_rodent, Portal_plots, by= c("plot", "month", "year"))

target_sp=c("PB", "PP", "PE", "PF", "PM", "RM")
target_repro=c("S", "M", "R")
target_month=c(1:12)


PB=portal1%>%
  filter(species=="PB",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

PP=portal1%>%
  filter(species=="PP",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

PE=portal1%>%
  filter(species=="PE",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

PF=portal1%>%
  filter(species=="PF",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

PM=portal1%>%
  filter(species=="PM",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

RM=portal1%>%
  filter(species=="RM",testes %in% target_repro, !(age=="J"), !is.na(wgt))%>%
  arrange(wgt)

total=rbind(PB,PP,PE,PF,RM,PM)

#only include adults and with weights > the minimum reproductive individual
target_ex=total%>%
  filter(treatment== "exclosure")%>%
  group_by(month, year)%>%
  summarize(count=n())

#all individuals counted per month

