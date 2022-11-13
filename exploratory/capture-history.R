require(reshape2)
#capture history
male1=portal_male%>%
  select(month, year, species, tag, testes, treatment)%>%
  reshape2::melt(id=c("tag", "treatment", "year","month", "species"))
dcast(tag+species+treatment+year~month, value.var = "variable", length)%>%
  arrange(year)

female1=portal_female%>%
  select(month, year, species, tag, nipples, pregnant, vagina, lactation, treatment)%>%
  reshape2::melt(id=c("tag", "treatment", "year","month", "species"))%>%
  mutate(reproductive=ifelse(value %in%c("P", "B", "S", "R", "E","L"), 1, 0))%>%
  select(-variable, -value)
  
female2=pivot_wider(female1, names_from=month, values_from=reproductive)

group_by(tag, month, year, species, treatment)%>%
  summarise(total_reprod=sum(reproductive))

  
  
dcast(tag+species+treatment+year~month, value.var = "reproductive", length)

 group_by(tag, species, treatment, year, month)%>%
  summarise(abundance=n(), total_reprod=sum(reproductive))
  dcast(tag+species+treatment+year~month, value.var = "reproductive", length)

  require(tidyr)

  female2=spread(female1, )  