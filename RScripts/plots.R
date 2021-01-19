all_fem_props=read.csv("./reconfigured_data/all_fem_props.csv")
all_props=read.csv("./reconfigured_data/all_props_NA.csv")
sum(all_fem_props$count)
sum(all_fem_props$repro, na.rm=T) #5234/15568=0.34

con1=all_fem_props%>%
  filter(trt=="control")
sum(ex1$repro, na.rm=T) #35%
sum(con1$count, na.rm=T)#31%

require(ggplot2)
require(ggplotgui)
ggplot(allmc, aes(x = year, y=proportion)) +
  geom_point(alpha=0.2)+
  geom_smooth(se = FALSE, method = 'lm') +
  #facet_grid( . ~ trt ) +
  theme_bw()+ggtitle("male (control)")


PEf=ggplot(PEf_props, aes(x = year, y = proportion, colour = trt)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  # facet_grid( . ~ month ) +
  theme_bw()+ggtitle("PE female")

PFf=ggplot(PFf_props, aes(x = year, y = proportion, colour = trt)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  # facet_grid( . ~ month ) +
  theme_bw()+ggtitle("PF female")

PMf=ggplot(PMf_props, aes(x = year, y = proportion, colour = trt)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  # facet_grid( . ~ month ) +
  theme_bw()+ggtitle("PM female")

PPf=ggplot(PPf_props, aes(x = year, y = proportion, colour = trt)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  # facet_grid( . ~ month ) +
  theme_bw()+ggtitle("PP female")

RMf=ggplot(RMf_props, aes(x = year, y = proportion, colour = trt)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  # facet_grid( . ~ month ) +
  theme_bw()+ggtitle("RM female")




PB_props%>%
  ggplot(mapping=aes(y=proportion, fill=trt))+
  geom_density(position = "dodge")+
  facet_wrap(~month)

# plots per treatment####
male_ex=all_sp%>%filter(!is.na(count), !is.na(repro))%>%group_by(Month,species, trt)%>%summarize(total_n=sum(count), total_rep=sum(repro), total_prop=total_rep/total_n)%>%select(Month, species, total_n, total_rep, total_prop, trt)%>%
  mutate(month=recode(Month, "Jan"="1", "Feb"="2", "Mar"="3", "Apr"="4",
     "May"="5","Jun"="6", "Jul"="7", "Aug"="8","Sept"= "9",
      "Oct"="10","Nov"="11", "Dec"="12"))


graph <- ggplot(male_ex, mapping=aes(x = reorder(Month,month), y = total_prop, colour = trt)) +
  geom_boxplot(notch = FALSE) +
theme_bw()
 graph

write.csv(male_ex, "monthly_male_trt.csv")


female_ex=all_female_sp%>%filter(!is.na(count), !is.na(repro))%>%group_by(Month,species, trt)%>%summarize(total_n=sum(count), total_rep=sum(repro), total_prop=total_rep/total_n)%>%select(Month, species, total_n, total_rep, total_prop, trt)%>%
  mutate(month=recode(Month, "Jan"="1", "Feb"="2", "Mar"="3", "Apr"="4",
                      "May"="5","Jun"="6", "Jul"="7", "Aug"="8","Sept"= "9",
                      "Oct"="10","Nov"="11", "Dec"="12"))

write.csv(female_ex, "monthly_female_trt.csv")

graph <- ggplot(female_ex, mapping=aes(x = reorder(Month,month), y = total_prop, colour = trt)) +
  geom_boxplot(notch = FALSE) +
  theme_bw()+ggtitle("Female")
graph
