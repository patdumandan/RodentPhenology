all_fem_props=read.csv("./reconfigured_data/all_fem_props.csv")

sum(all_fem_props$count)
sum(all_fem_props$repro, na.rm=T) #5234/15568=0.34

con1=all_fem_props%>%
  filter(trt=="control")
sum(ex1$repro, na.rm=T) #35%
sum(con1$count, na.rm=T)#31%

require(ggplot2)
require(ggplotgui)
ggplot(all_props, aes(x = year, y = proportion, colour = Month)) +
  geom_point()+
  geom_smooth(se = FALSE, method = 'lm') +
  #facet_grid( . ~ month ) +
  theme_bw()+ggtitle("PP male")


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

qplot(y=proportion, data=PBf_props, geom="histogram")
