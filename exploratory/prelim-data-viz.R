#preliminary visualization 
#try after running data_curation_final.R
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
