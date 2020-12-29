library(dplyr)
library(tidyr)
Portal_rodent=read.csv("./PortalData/Rodents/Portal_rodent.csv")
Portal_plots=read.csv("./PortalData/SiteandMethods/Portal_plots.csv")
portal1=inner_join(Portal_rodent, Portal_plots, by= c("plot", "month", "year"))

target_sp=c("PB", "PP", "PE", "PF", "PM", "RM")
target_repro=c("S", "M", "R")
target_month=c(1:12)


PB=portal1%>%
  filter(species=="PB",testes %in% target_repro,wgt>=16)%>%
  arrange(wgt)

PP=portal1%>%
  filter(species=="PP",testes %in% target_repro, wgt>=10)%>%
  arrange(wgt)

PE=portal1%>%
  filter(species=="PE",testes %in% target_repro, wgt >=7)%>%
  arrange(wgt)

PF=portal1%>%
  filter(species=="PF",testes %in% target_repro, wgt>=4)%>%
  arrange(wgt)

PM=portal1%>%
  filter(species=="PM",testes %in% target_repro, wgt >=11)%>%
  arrange(wgt)
#1 individual identified as juv but size within reproductive??

RM=portal1%>%
  filter(species=="RM",testes %in% target_repro, wgt>=4)%>%
  arrange(wgt)

total_repro=rbind(PB,PP,PE,PF,RM,PM)

####ALL TARGET SPECIES ACROSS ALL YEARS####
#only include adults and with weight > the minimum reproductive individual
target_ex=total_repro%>%
  filter(treatment== "exclosure")%>%
  group_by(month)%>%
  summarize(reproductive_ex=n())

target_con=total_repro%>%
  filter(treatment== "control")%>%
  group_by(month)%>%
  summarize(reproductive_con=n())

#all individuals counted per month
total_ex=portal1%>%
      filter(species %in% target_sp, treatment=="exclosure")%>%
      group_by(month)%>%
      summarize(total=n())

total_con=portal1%>%
  filter(species %in% target_sp, treatment=="control")%>%
  group_by(month)%>%
  summarize(total_cont=n())

prop_ex=inner_join(target_ex, total_ex)%>%
  mutate(proportion=reproductive_ex/total*100, treatment="exclosure")

prop_con=inner_join(target_con, total_con)%>%
  mutate(proportion_con=reproductive_con/total_cont*100, treatment="control")

prop_con$mon=prop_con$month

prop_ex%>%
  ggplot(mapping=aes(x=month, y=proportion))+
  geom_col(position=position_dodge())+
  scale_x_continuous(breaks=c(1:12),
    labels=c("1"="Jan", "2"="Feb", "3"="Mar",
                            "4"="Apr", "5"= "May", "6"= "Jun",
                            "7"="Jul", "8"="Aug", "9"= "Sept",
                            "10"="Oct", "11"="Nov", "12"="Dec"))+
  xlab("Month")+ylab("Proportion of reproductive individuals")

prop_con%>%
  ggplot(mapping=aes(x=month, y=proportion_con))+
  geom_col(position=position_dodge())+
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1"="Jan", "2"="Feb", "3"="Mar",
                              "4"="Apr", "5"= "May", "6"= "Jun",
                              "7"="Jul", "8"="Aug", "9"= "Sept",
                              "10"="Oct", "11"="Nov", "12"="Dec"))+
  xlab("Month")+ylab("Proportion of reproductive individuals")

#### PB proportions#####

PBex=portal1%>%
  filter(species=="PB",testes %in% target_repro,wgt>=16, treatment=="exclosure")%>%
  group_by(month, year)%>%
  summarize(repro=n())%>%
  mutate(trt="exclosure")

PBcon=portal1%>%
  filter(species=="PB",testes %in% target_repro,wgt>=16, treatment=="control")%>%
  group_by(month, year)%>%
  summarize(repro=n())%>%
  mutate(trt="control")

PB_all_ex=portal1%>%
  filter(species=="PB", treatment=="exclosure", sex=="M")%>%
  group_by(month, year)%>%
  summarize(count=n())%>%
  mutate(trt="exclosure")

PB_all_con=portal1%>%
  filter(species=="PB", treatment=="control", sex=="M")%>%
  group_by(month, year)%>%
  summarize(count=n())%>%
  mutate(trt="control")

df1=left_join(PB_all_ex, PBex)
df2=left_join(PB_all_con,PBcon)
df3=rbind(df1,df2)

df4=df3%>%
  mutate(proportion=repro/count)
write.csv(df4, "PB_props.csv")
PB_props=read.csv("PB_props.csv")

DF=df%>%
  mutate(proportion=repro/count*100)
