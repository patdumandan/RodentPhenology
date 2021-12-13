#these are plots for the derivatives###
#to plot estimated smooths, remove the lines of code calculating derivatives (e.g., p1d=derivatives(model))
#PBs####

tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf)
source(tmpf)

#female####

#control####
p1=ggplot(PB_female_con, aes(y=proportion, x=month, col=treatment)) +
  geom_point() +
  ylab("P(breeding)")+
  #geom_smooth(method = 'gam', formula = y ~ s(x, bs="cc"))+
  ggtitle("Bailey's pocket mouse")+ facet_wrap(~sex)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p1d=derivatives(pbf_con3)
p1=draw(p1d, select=1, shade=F, jit=T)

p1=draw(pbf_con3, select=1, shade=F, jit=T) #for estimated smooth

p1=p1+annotate("rect", xmin=1, xmax=3, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=6, xmax=8, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("A")
 
p1=p1%>% annotate_figure(left=text_grob("female",face = "bold"), top=text_grob("control",face = "bold"))
                         

#exclosure####
p2d=derivatives(pbf_ex3)

p2=draw(p2d, select=1, shade=F, jit=T)

p2=draw(pbf_ex3, select=1, shade=F, jit=T) #for estimated smooth

p2=p2+annotate("rect", xmin=2, xmax=4, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=7, xmax=9, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("B")

p2=p2%>% annotate_figure(top=text_grob("exclosure", face = "bold"))

#PBF=ggarrange(p1,p2)%>%annotate_figure(top=text_grob("Bailey's pocket mouse"))

#male####
#control####
p3d=derivatives(pbm_con3)

p3=draw(p3d, select=1, shade=F, jit=T)

p3=draw(pbm_con3, select=1, shade=F, jit=T) #for estimated smooth

p3=p3+scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("C")

p3=p3%>% annotate_figure(left=text_grob("male", face = "bold"))

#exclosure####
p4d=derivatives(pbm_ex3)

p4=draw(p4d, select=1, shade=F, jit=T)

p4=draw(pbm_ex3, select=1, shade=F, jit=T) #for estimated smooth

p4=p4+annotate("rect", xmin=1, xmax=3, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=7, xmax=7.5, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("D")


PB=ggarrange(p1,p2,p3, p4, ncol=2, nrow=2, common.legend = TRUE)
  

#PPs####

#female####

#control####
p5d=derivatives(ppf_con3)

p5=draw(p5d, select=1, shade=F, jit=T)

p5=draw(ppf_con3, select=1, shade=F, jit=T) #for estimated smooth

p5=p5+annotate("rect", xmin=11, xmax=12, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=1, xmax=1.5, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=5, xmax=8, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("A")

p5=p5%>% annotate_figure(left=text_grob("female",face = "bold"), top=text_grob("control",face = "bold"))

#exclosure####
p6d=derivatives(ppf_ex3)

p6=draw(p6d, select=1, shade=F, jit=T)

p6=draw(ppf_ex3, select=1, shade=F, jit=T) #for estimated smooth


p6=p6+annotate("rect", xmin=10, xmax=12, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=5, xmax=8, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("B")

p6=p6%>% annotate_figure(top=text_grob("exclosure", face = "bold"))

#male####

#control####
p7d=derivatives(ppm_con3)

p7=draw(p7d, select=1, shade=F, jit=T)

p7=draw(ppm_con3, select=1, shade=F, jit=T) #for estimated smooth

p7=p7+annotate("rect", xmin=1, xmax=3, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=6, xmax=8, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("C")

p7=p7%>% annotate_figure(left=text_grob("male", face = "bold"))

#exclosure####
p8d=derivatives(ppm_ex3)

p8=draw(p8d, select=1, shade=F, jit=T)

p8=draw(ppm_ex3, select=1, shade=F, jit=T) #for estimated smooth

p8=p8+annotate("rect", xmin=1, xmax=2, ymin=-Inf, ymax=Inf,fill="blue", alpha=0.1)+
  annotate("rect", xmin=7, xmax=8, ymin=-Inf, ymax=Inf,fill="red", alpha=0.1)+
  scale_x_discrete(name="month", limits=c("Jan", "Feb", "Mar", "Apr","May", "Jun",
                                          "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme_classic()+ggtitle("D")


PP=ggarrange(p5,p6,p7, p8, ncol=2, nrow=2, common.legend = TRUE)
