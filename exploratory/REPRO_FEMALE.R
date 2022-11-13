require(viridis)
require(portalr)
require(dplyr)
require(ggplot2)


#FULL DATASET####

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
  quiet = FALSE)%>%
  filter(!is.na(sex)& !is.na(treatment))

#SELECT FEMALES IN CONTROL AND EXCLOSURE ONLY####

portal_fem=Portal_data_indiv%>%
  filter(treatment==c("exclosure", "control")&
           sex=="F")

#13854 records of F in exclosure and control

#SELECT REPRODUCTIVE FEMALES ONLY####

repro_female=portal_fem%>%
  filter(vagina==c("S", "P", "B")| pregnant=="P" | nipples==c("R", "E", "B") | lactation=="L")

#17% or 2341 REPRODUCTIVE OF 14031 INDIVS

repro_female=as.data.frame(repro_female)
write.csv(repro_female, "repro_female.csv")

#PB-SPECIFIC####

#SIZE THRESHOLDS####
PBf_tot=portal_fem%>%filter(species=="PB")

#2598 TOTAL PB females

PBf_repro=repro_female%>%
  filter(species=="PB" & wgt >=23)%>%
  group_by(month, year, treatment)%>%
  summarise(reproductive=n()) 

#295 reproductive or 11%

#CALCULATE PROPORTIONS####

#exclosure####
PBf_repro_ex=PBf_repro%>%filter(treatment=="exclosure") 

#217 of 295 or 74% breeding F in exclosure

PBf_all_ex=PBf_tot%>%
  filter(treatment=="exclosure")%>%
  group_by(month,year)%>%
  summarize(count=n()) 

#1850 OF 2598 PBs or 71% in exclosure 

#PBfex=portal_fem%>%filter(species=="PB"&treatment=="exclosure")

PBf_total_ex=left_join(PBf_all_ex, PBf_repro_ex)%>%
  mutate(reproductive=replace_na(reproductive, 0), treatment=replace_na(treatment, "exclosure"),
         proportion=reproductive/count)%>%
  arrange(proportion)


#control####
PBf_repro_con=PBf_repro%>%filter(treatment=="control") 

#78 of 295 or 26% breeding F in control

PBf_all_con=PBf_tot%>%
  filter(treatment=="control")%>%
  group_by(month,year)%>%
  summarize(count=n()) 

#748 OF 2598 PBs or 29% in control 

#PBfcon=portal_fem%>%filter(species=="PB",treatment=="control")

PBf_total_con=left_join(PBf_all_con, PBf_repro_con)%>%
  mutate(reproductive=replace_na(reproductive, 0), treatment=replace_na(treatment, "control"),
         proportion=reproductive/count)%>%
  arrange(proportion)

Pbf_dat=rbind(PBf_total_con, PBf_total_ex)
Pbf_dat=as.data.frame(Pbf_dat)
write.csv(Pbf_dat, "Pb_females.csv")

ggplot(reprod_PB_male, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  stat_smooth(method = 'gam', formula = y ~ s(x))+
  ggtitle("PB females")+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                        "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
