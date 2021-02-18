tot_dat=read.csv("./reconfigured_data/raw_cleaned/reprod_propn_male.csv")
require(dplyr)
PB_dat_M=tot_dat%>%filter(species=="PB", !(treatment=="spectabs"))
PB_dat_M[is.na(PB_dat_M)] <- 0 #set non-detects to 0
PB_dat_M$trt<-ifelse(PB_dat_M$treatment=="control", 0, 1) 
PB_dat_M$years=(PB_dat_M$year-mean(PB_dat_M$year))/(2*sd(PB_dat_M$year)) #standardize year
PBprop=PB_dat_M$proportion
PB_dat_M$mon_cos= cos(2*pi*(PB_dat_M$month/12))
PB_dat_M$mon_sin= sin(2*pi*(PB_dat_M$month/12))
PB_dat_M$spcode=as.integer(PB_dat_M$species)


dat_list=list(
  N=length(PB_dat_M$month),
  y=PB_dat_M$reproductive,
  n=PB_dat_M$abundance,
  year=PB_dat_M$years,
  treatment=PB_dat_M$trt,
  month=as.integer(PB_dat_M$month),
  mon_cos=PB_dat_M$mon_cos, 
  mon_sin=PB_dat_M$mon_sin,
  Nmon=length(unique(PB_dat_M$month)),
  Nsp=length(unique(PB_dat_M$species)))


y=PB_dat_M$proportion
yrep2=rstan::extract(mod_int9)$pred_y
con_pb=yrep2[,which(PB_dat_M$treatment=="control"& PB_dat_M$month==3)]
con_pbmat=as.matrix(con_pb)
con_pbs=con_pbmat[1:300,]
matplot(t(con_pbs), type="l", col="grey", main="PB control (March)")
mean_con_pb=apply(con_pb, 2, mean)
con_pb_obs=PB_dat_M%>%filter(treatment=="control"& month==3)
lines(mean_con_pb~c(1:length(mean_con_pb)), col="white")
points(con_pb_obs$proportion, col="black", cex=1 )