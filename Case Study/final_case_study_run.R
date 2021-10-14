

library(Rstrap)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(RColorBrewer)
library(rasterVis)

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Data\\DFO")

## Load survey data compiled for Rstrap
load("Spatial_RV_Data/converted_set_details_2019-04-23.Rdata")
load("Spatial_RV_Data/converted_length-frequency_data_2019-04-23.Rdata")
load("Spatial_RV_Data/set_details_2019-04-23.Rdata")
load("Spatial_RV_Data/length-frequency_data_2019-04-23.Rdata")
load("Spatial_RV_Data/age-growth_data_2019-04-23.Rdata")

## Subset data to sand lance
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 694)), ]
con.lf <- con.lf[con.lf$spec == 694, ]
setdet <- setdet[(setdet$rec == 5 | (setdet$rec == 6 & setdet$spec == 694)), ]
lf <- lf[lf$spec == 694, ]
ag <- ag[ag$spec == 694, ]

## Convert prey codes to text
prey_codes <- read.csv("Diet_specs/stomach_contents_codes.csv", stringsAsFactors = FALSE)
names(prey_codes) <- c("code", "name", "group")
ind <- !(is.na(ag$stomach) | nchar(ag$stomach) != 5 |
           grepl("[[:alpha:]|[:blank:]|[:punct:]]", ag$stomach))
ag$fullness[ind] <- as.numeric(substr(ag$stomach[ind], 1, 1)) / 9
ag$empty <- ag$fullness == 0
ag$prey1[ind] <- as.numeric(substr(ag$stomach[ind], 2, 3))
ag$prey2[ind] <- as.numeric(substr(ag$stomach[ind], 4, 5))
ag$prey1 <- prey_codes$name[match(ag$prey1, prey_codes$code)]
ag$prey2 <- prey_codes$name[match(ag$prey2, prey_codes$code)]
ag$prey1[ag$empty] <- "Empty"
ag$prey2[ag$empty] <- "Empty"
ag$prey1_group <- prey_codes$group[match(ag$prey1, prey_codes$name)]
ag$prey2_group <- prey_codes$group[match(ag$prey2, prey_codes$name)]

## Save for raw data explorations
rv_data <- list(setdet = setdet)
max_year <- max(setdet$survey.year)


# Run DFO program RStrap to access all cod data in 3LNO from 1995 onwards
fall_camp<- strat.fun(setdet = rv_data$setdet,
                      program = "strat2",
                      data.series = c("Campelen"),
                      species = 694,
                      survey.year = 1995:max_year,
                      season = "spring",
                      NAFOdiv = c("3L","3N","3O"),
                      sex = c("male", "female", "unsexed"),
                      plot.results = FALSE,
                      export = NULL)

#extract the set details
fall_camp_df<-fall_camp$raw.data$set.details

#change number caught to binary
fall_camp_df$pos_catch <- ifelse(fall_camp_df$number>0, 1, 0)

#change name of dataframe
just_onespp<-fall_camp_df

#extract the years in that dataframe
uni_yrs<-unique(just_onespp$year)

#create column with unique set ids
just_onespp$map<-paste(just_onespp$vessel, just_onespp$trip, just_onespp$set, sep=".")

#create spatial dataframe of positive catches
anomaly_sp<-data.frame(long=-just_onespp$long.start, lat=just_onespp$lat.start)
anomaly_spdf<-SpatialPointsDataFrame(anomaly_sp, data.frame(pa=just_onespp$pos_catch))
anomaly_spdf2<-SpatialPointsDataFrame(anomaly_sp, data.frame(pa=just_onespp$pos_catch, map=just_onespp$map))

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Data\\Spatial Data")

#import spatial data about the grand bank
nafo <- readOGR(dsn = "Divisions", layer = "Divisions", verbose=FALSE)
land <- readOGR(dsn = "ne_10m_bathymetry_L_0", layer = "ne_10m_bathymetry_L_0", verbose=FALSE)
shelf <- readOGR(dsn = "ne_10m_bathymetry_K_200", layer = "ne_10m_bathymetry_K_200", verbose=FALSE)

# now I define a raster based on the extent of the Grand Bank:
r = raster(extent(-55, -46.3, 42, 49.25))

res(r)=c(0.5,0.5) #resolution of the raster (0.5 x 0.5 deg)

# crop shelf polygon data to the same extent as the raster
shelf2<-crop(shelf,(extent(c(-55,-46.3,42.5,49.15))))

#make sure the CRS's match
crs(anomaly_spdf)<-crs(shelf2)
crs(anomaly_spdf2)<-crs(shelf2)

#rasterize the mean number of positive sand lance catches across the grand bank
anomaly_rast<-rasterize(anomaly_spdf,r, fun=mean, na.rm=TRUE)
anomaly_rast<-dropLayer(anomaly_rast,1) 

crs(anomaly_rast)<-crs(anomaly_spdf)

# define plotting parameters
myTheme <- BTCTheme()
myTheme$panel.background$col = 'lightgrey' 
p.strip <- list(cex=1.5, lines=1)

#change areas with a mean catch of 0 to NA to remove samples from those areas later on
#could modify this threshold to remove samples at different levels
anomaly_rast[anomaly_rast==0] = NA

#create a spatialpolygons object from the raster grid for plotting
grid_layer <- as(r,'SpatialPolygonsDataFrame')

#create polygons object from the raster
rastpoly<-rasterToPolygons(anomaly_rast)

#remove na's
rastpoly@data <- na.omit(rastpoly@data)

#intersect locations that aren't NA with the spdf. Essentially removes all locations that had samples where mean catch was 0
new_thing<-raster::intersect(anomaly_spdf2, rastpoly)

library(TMB)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Dissertation Plans\\Chapter 2\\Final Products\\NLFPM\\Models")

#prep the code to run the different models
dyn.unload("LFPM")
compile("LFPM.cpp")
dyn.load("LFPM")

dyn.unload("NLFPM")
compile("NLFPM.cpp")
dyn.load("NLFPM")

dyn.unload("both_species_model")
compile("both_species_model.cpp")
dyn.load("both_species_model")

#load all data
load("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Dissertation Plans\\Chapter 2\\Stomach-Content-Simulation\\Noels_parm\\ch2_data_7_26.RData")

#subset the data for species of interest
ampl_call<-subset(ch2_data$spring_call_dat, name=="american plaice")
ampl_full<-subset(ch2_data$spring_sto_dat, Predator_spp_name=="plaice")

cod_call<-subset(ch2_data$spring_call_dat, name=="cod,atlantic")
cod_full<-subset(ch2_data$spring_sto_dat, Predator_spp_name=="cod")

spring_camp_surv<-ch2_data$spring_camp_surv

#subset out samples in locations that dont catch sand lance
spring_camp_surv<-spring_camp_surv[spring_camp_surv$map %in% new_thing@data$map,]

ampl_call<-ampl_call[ampl_call$map %in% new_thing@data$map,]
ampl_full<-ampl_full[ampl_full$map %in% new_thing@data$map,]
cod_call<-cod_call[cod_call$map %in% new_thing@data$map,]
cod_full<-cod_full[cod_full$map %in% new_thing@data$map,]

#subset out bad years and to account for ontogenetic diet
ampl_call<-subset(ampl_call,length>25 & survey.year!=2006 & survey.year!=2016 & survey.year!=2017 & survey.year!=2018)
ampl_sto<-subset(ampl_full, LENGTH>25 & Year>2012)

cod_call<-subset(cod_call, length>25)
cod_sto<-subset(cod_full, LENGTH>25 & Year<1998|Year>2012)

length(spring_camp_surv$survey.year)
length(ampl_call$year)
length(ampl_full$ID_pred)
length(cod_call$year)
length(cod_full$ID_pred)


#prep data for model
new_df<-data.frame(pa=c(spring_camp_surv$trawl_pa, ampl_call$pa, ampl_sto$pa),
                   year=c(spring_camp_surv$survey.year, ampl_call$survey.year, ampl_sto$Year),
                   idex=c(rep(0, length(spring_camp_surv$survey.year)), rep(1, length(ampl_call$survey.year)),
                          rep(2, length(ampl_sto$Year))))

#double check that bad years arent there
new_df<-subset(new_df, year!=1995 & year!=2019)

#create objects to compare model output to
mean_trawl<-aggregate(spring_camp_surv$trawl_pa, by=list(spring_camp_surv$year), FUN=mean, na.rm=TRUE)
mean_call<-aggregate(ampl_call$pa, by=list(ampl_call$year), FUN=mean, na.rm=TRUE)
mean_sto<-aggregate(ampl_sto$pa, by=list(ampl_sto$Year), FUN=mean, na.rm=TRUE)
colnames(mean_call)<-c("Group.1", "biased")
colnames(mean_sto)<-c("Group.1", "biased")

year<-sort(unique(new_df$year))

###############################################3
#######################################
############American plaice models
#######################################
#############################################

#prep data list
tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=3,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  pa=new_df$pa)

#prep parameters
parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lk=0
)


###RUNNING A TYPE 1 MODEL
obj <- MakeADFun(tmb.data,parameters,DLL="LFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)

#save output
t1_ampl_rep = obj$report()

t1_ampl_aic<-2*opt$objective+ 2*length(opt$par)#aic

#prep data and parms for NLFPM
tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=3,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  k=1,
  pa=new_df$pa
)

parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lbeta = log(2),
  lchi = log(0.5),
  logrw_var=log(1)
)

# Run NLFPM
obj <- MakeADFun(tmb.data,parameters,DLL="NLFPM",inner.control=list(maxit=2000,trace=F),
                 random=c("iye"),
                 silent=TRUE)
opt<-nlminb(obj$par,obj$fn,obj$gr, control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)

t2_ampl_aic<-2*opt$objective+ 2*length(opt$par)

#save output
ampl_rep<-obj$report() 
ampl_sdrep<-sdreport(obj) 
save(ampl_sdrep, file="ampl_sdrep_nozeroes.RData")

#plotting objects from output
mean_call$linear<-((mean_call$biased*(ampl_rep$chi^ampl_rep$beta))/(1-mean_call$biased))^(1/ampl_rep$beta)
mean_sto$linear<-((mean_sto$biased*(ampl_rep$chi^ampl_rep$beta))/(1-mean_sto$biased))^(1/ampl_rep$beta)

joined_df<-plyr::join(mean_trawl, mean_call, by="Group.1")
ampl_joined_df2<-plyr::join(joined_df, mean_sto, by="Group.1")
colnames(ampl_joined_df2)<-c("Year","Trawl","Raw_call","Fix_call","Raw_sto","Fix_sto")

ampl_joined_df2[ampl_joined_df2<0.01] <- NA

newx<-seq(from=0, to=1, by=0.001)
K<-1
chi<-ampl_rep$chi
beta<-ampl_rep$beta
rate<-((K*(newx^beta))/((chi^beta)+(newx^beta)))


chi_high=ampl_rep$chi+qnorm(0.975)*ampl_sdrep$sd[25]
chi_low=ampl_rep$chi-qnorm(0.975)*ampl_sdrep$sd[25]

beta_high=ampl_rep$beta+qnorm(0.975)*ampl_sdrep$sd[24]
beta_low=ampl_rep$beta-qnorm(0.975)*ampl_sdrep$sd[24]

rate_high<-(newx^beta_high)/(chi_high^beta_high + newx^beta_high)
rate_low<-(newx^beta_low)/(chi_low^beta_low + newx^beta_low)

curve_df<-data.frame(x=newx, rate=rate, rate_low=rate_low, rate_high=rate_high)

ampl_curve<-ggplot()+
  geom_ribbon(data=curve_df, aes(x=x, y=rate, ymin=rate_low, ymax=rate_high), alpha=0.1)+
  geom_line(data=curve_df, aes(x=x, y=rate), size=2)+
  geom_point(data=ampl_joined_df2, aes(x=Trawl, y=Raw_call, col="#92c5de"), size=3, alpha=0.7)+
  geom_point(data=ampl_joined_df2, aes(x=Trawl, y=Raw_sto, col="#f4a582"), size=3, alpha=0.7)+
  geom_segment(aes(x=0.7,y=0.7, xend=0.75, yend=0.52),
               lineend = "round", linejoin = "round",
               size = 1.5, arrow = arrow(length = unit(0.2, "inches")))+
  scale_color_manual(name="Data",labels=c("Called Stomach Contents","Full Stomach Contents"),values=c("#92c5de", "#f4a582"),
                     guide = guide_legend(override.aes=aes(fill=NA)))+
  ggtitle("American plaice")+
  xlim(0,1)+ylim(0,1)+
  xlab("")+ylab("")+
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right",
        legend.text=element_text(size=16),legend.key = element_rect(fill = NA),legend.title=element_text(size=16, face="bold"),
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))+
  annotate(geom="text",x=0.7, y=0.8, label=expression( frac(italic(pt[y])^2.98, 0.78^2.98~+italic(pt[y])^2.98)), size=6)


low_rho = ampl_sdrep$value - qnorm(0.975)*ampl_sdrep$sd
high_rho = ampl_sdrep$value + qnorm(0.975)*ampl_sdrep$sd

ts_df<-data.frame(ts=ampl_rep$iye, time=year, low_sd=low_rho, high_sd=high_rho, call=log(1/(1-ampl_joined_df2$Fix_call)),
                  trawl=log(1/(1-ampl_joined_df2$Trawl)), sto=log(1/(1-ampl_joined_df2$Fix_sto)))

ampl_ts<-ggplot()+
  geom_ribbon(data=ts_df, aes(x=year, y=ts, ymin=low_sd, ymax=high_sd), alpha=0.2)+
  geom_line(data=ts_df, aes(x=time, y=ts, colour="1"),size=1.5)+
  geom_line(data=ts_df, aes(x=time, y=call,colour="2"), linetype="dashed", size=1.4)+
  geom_line(data=ts_df, aes(x=time, y=sto, colour="3"), linetype="dashed", size=1.4)+
  geom_line(data=ts_df, aes(x=time, y=trawl, colour="4"), linetype="dashed", size=1.4)+
  scale_color_manual(name="Single Species Models",values=c("#656565","#035d8a","#f4a582","#c2a5cf"),
                     labels=c("Estimated N","Called Stomach","Full Stomach","Trawl"))+
  ggtitle("American plaice")+
  xlab("")+ylab("")+
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right",
        legend.text=element_text(size=16),legend.key = element_rect(fill = NA),legend.title=element_text(size=16, face="bold"),
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))


###############################################3
#######################################
############Atlantic cod models
#######################################
#############################################

new_df<-data.frame(pa=c(spring_camp_surv$trawl_pa, cod_call$pa, cod_sto$pa),
                   year=c(spring_camp_surv$survey.year, cod_call$survey.year, cod_sto$Year),
                   idex=c(rep(0, length(spring_camp_surv$survey.year)), rep(1, length(cod_call$survey.year)),
                          rep(2, length(cod_sto$Year))))

new_df<-subset(new_df, year!=1995 & year!=2019)

mean_trawl<-aggregate(spring_camp_surv$trawl_pa, by=list(spring_camp_surv$year), FUN=mean, na.rm=TRUE)
mean_call<-aggregate(cod_call$pa, by=list(cod_call$year), FUN=mean, na.rm=TRUE)
mean_sto<-aggregate(cod_sto$pa, by=list(cod_sto$Year), FUN=mean, na.rm=TRUE)
colnames(mean_call)<-c("Group.1", "biased")
colnames(mean_sto)<-c("Group.1", "biased")

year<-sort(unique(new_df$year))

tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=3,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  pa=new_df$pa)

parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lk=0
)


###RUNNING A TYPE 1 MODEL
obj <- MakeADFun(tmb.data,parameters,DLL="LFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
t1_cod_rep = obj$report()

t1_cod_aic<-2*opt$objective+ 2*length(opt$par)#aic


tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=3,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  k=1,
  pa=new_df$pa
)

parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lbeta = log(2),
  lchi = log(1),
  logrw_var=log(1)
)


obj <- MakeADFun(tmb.data,parameters,DLL="NLFPM",inner.control=list(maxit=2000,trace=F),
                 random=c("iye"),
                 silent=TRUE)

opt<-nlminb(obj$par,obj$fn,obj$gr, control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)

obj$gr()

t2_cod_aic<-2*opt$objective+ 2*length(opt$par)

cod_rep<-obj$report()
cod_sdrep<-sdreport(obj)

mean_call$linear<-((mean_call$biased*(cod_rep$chi^cod_rep$beta))/(1-mean_call$biased))^(1/cod_rep$beta)
mean_sto$linear<-((mean_sto$biased*(cod_rep$chi^cod_rep$beta))/(1-mean_sto$biased))^(1/cod_rep$beta)

joined_df<-plyr::join(mean_trawl, mean_call, by="Group.1")
cod_joined_df2<-plyr::join(joined_df, mean_sto, by="Group.1")
colnames(cod_joined_df2)<-c("Year","Trawl","Raw_call","Fix_call","Raw_sto","Fix_sto")

cod_joined_df2[cod_joined_df2<0.01] <- NA


newx<-seq(from=0, to=1, by=0.001)
K<-1
chi<-cod_rep$chi
beta<-cod_rep$beta
rate<-((K*(newx^beta))/((chi^beta)+(newx^beta)))

chi_high=cod_rep$chi+qnorm(0.975)*cod_sdrep$sd[25]
chi_low=cod_rep$chi-qnorm(0.975)*cod_sdrep$sd[25]

beta_high=cod_rep$beta+qnorm(0.975)*cod_sdrep$sd[24]
beta_low=cod_rep$beta-qnorm(0.975)*cod_sdrep$sd[24]

rate_high<-(newx^beta_high)/(chi_high^beta_high + newx^beta_high)
rate_low<-(newx^beta_low)/(chi_low^beta_low + newx^beta_low)

curve_df<-data.frame(x=newx, rate=rate, rate_low=rate_low, rate_high=rate_high)

cod_curve<-ggplot()+
  geom_ribbon(data=curve_df, aes(x=x, y=rate, ymin=rate_low, ymax=rate_high), alpha=0.1)+
  geom_line(data=curve_df, aes(x=x, y=rate), size=2)+
  geom_point(data=cod_joined_df2, aes(x=Trawl, y=Raw_call), col="#92c5de", size=3, alpha=0.7)+
  geom_point(data=cod_joined_df2, aes(x=Trawl, y=Raw_sto), col="#f4a582", size=3, alpha=0.7)+
  geom_segment(aes(x=0.55,y=0.75, xend=0.6, yend=0.54),
               lineend = "round", linejoin = "round",
               size = 1.5, arrow = arrow(length = unit(0.2, "inches")))+
  ggtitle("Atlantic cod")+
  xlim(0,1)+ylim(0,1)+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))+
  annotate(geom="text",x=0.55, y=0.85, label=expression( frac(italic(pt[y])^5.50, 0.62^5.50~+italic(pt[y])^5.50)), size=6)

low_rho = cod_sdrep$value - qnorm(0.975)*cod_sdrep$sd
high_rho = cod_sdrep$value + qnorm(0.975)*cod_sdrep$sd

ts_df<-data.frame(ts=cod_rep$iye, time=year, low_sd=low_rho, high_sd=high_rho, call=log(1/(1-cod_joined_df2$Fix_call)),
                  trawl=log(1/(1-cod_joined_df2$Trawl)), sto=log(1/(1-cod_joined_df2$Fix_sto)))

cod_ts<-ggplot()+
  geom_ribbon(data=ts_df, aes(x=year, y=ts, ymin=low_sd, ymax=high_sd), alpha=0.2)+
  geom_line(data=ts_df, aes(x=time, y=ts, colour="1"),size=1.5)+
  geom_line(data=ts_df, aes(x=time, y=call,colour="2"), linetype="dashed", size=1.4)+
  geom_line(data=ts_df, aes(x=time, y=sto, colour="3"), linetype="dashed", size=1.4)+
  geom_line(data=ts_df, aes(x=time, y=trawl, colour="4"), linetype="dashed", size=1.4)+
  scale_color_manual(name="",values=c("#656565","#035d8a","#f4a582","#c2a5cf"),
                     labels=c("Estimated N","Called Stomach","Full Stomach","Trawl"))+
  ggtitle("Atlantic cod")+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))

###############################################3
#######################################
############Both species models
#######################################
#############################################

new_df<-data.frame(pa=c(spring_camp_surv$trawl_pa, ampl_call$pa, ampl_sto$pa, cod_call$pa, cod_sto$pa),
                   year=c(spring_camp_surv$survey.year,ampl_call$survey.year, ampl_sto$Year, cod_call$survey.year, cod_sto$Year),
                   idex=c(rep(0, length(spring_camp_surv$survey.year)),rep(1, length(ampl_call$survey.year)),
                          rep(2, length(ampl_sto$Year)),rep(3, length(cod_call$survey.year)),
                          rep(4, length(cod_sto$Year))))

new_df<-subset(new_df, year!=1995 & year!=2019)

tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=5,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  pa=new_df$pa)

parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lk=0
)


###RUNNING A TYPE 1 MODEL
obj <- MakeADFun(tmb.data,parameters,DLL="LFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
t1_all_rep = obj$report()

t1_all_aic<-2*opt$objective+ 2*length(opt$par)#aic


tmb.data<-list(
  n=length(new_df$pa),
  nyrs=length(unique(new_df$year)),
  ndex=5,
  iyear=new_df$year-min(new_df$year),
  idex=new_df$idex,
  k=1,
  pa=new_df$pa
)

parameters<-list(
  iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
  lbeta = rep(log(2),2),
  lchi = rep(log(0.5),2)
)

obj <- MakeADFun(tmb.data,parameters,DLL="both_species_model",inner.control=list(maxit=2000,trace=F),
                 random=c("iye"),
                 silent=TRUE)

opt<-nlminb(obj$par,obj$fn,obj$gr, control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)

obj$gr()

t2_all_aic<-2*opt$objective+ 2*length(opt$par)

all_rep<-obj$report()
all_sdrep<-sdreport(obj)

newx<-seq(from=0, to=1, by=0.001)
K<-1
chi<-all_rep$chi[1]
beta<-all_rep$beta[1]
ampl_rate<-((K*(newx^beta))/((chi^beta)+(newx^beta)))
chi<-all_rep$chi[2]
beta<-all_rep$beta[2]
cod_rate<-((K*(newx^beta))/((chi^beta)+(newx^beta)))

curve_df<-data.frame(x=newx, ampl_rate=ampl_rate, cod_rate=cod_rate)

both_curve<-ggplot()+
  geom_line(data=curve_df, aes(x=x, y=ampl_rate), size=2, colour="#035d8a")+
  geom_line(data=curve_df, aes(x=x, y=cod_rate), size=2, colour="#f4a582",)+
  geom_segment(aes(x=0.57,y=0.72, xend=0.63, yend=0.52),
               lineend = "round", linejoin = "round",
               size = 1.5, arrow = arrow(length = unit(0.2, "inches")), colour="#f4a582")+
  geom_segment(aes(x=0.78,y=0.2, xend=0.7, yend=0.3),
               lineend = "round", linejoin = "round",
               size = 1.5, arrow = arrow(length = unit(0.2, "inches")), colour="#035d8a")+
  ggtitle("Both species")+
  xlim(0,1)+ylim(0,1)+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))+
  annotate(geom="text",x=0.75, y=0.1, label=expression(paste("Plaice =",frac(italic(pt[y])^2.20, 0.97^2.20~+italic(pt[y])^2.20))), size=6, colour="#035d8a")+
  annotate(geom="text",x=0.5, y=0.84, label=expression(paste("Cod =", frac(italic(pt[y])^4.62, 0.67^4.62~+italic(pt[y])^4.62))), size=6, colour="#f4a582")

low_rho = all_sdrep$value - qnorm(0.975)*all_sdrep$sd
high_rho = all_sdrep$value + qnorm(0.975)*all_sdrep$sd

ts_df<-data.frame(ts=all_rep$iye, time=year, low_sd=low_rho, high_sd=high_rho)

cod_ts_df<-data.frame(ts=cod_rep$iye, time=year)
ampl_ts_df<-data.frame(ts=ampl_rep$iye, time=year)


both_ts<-ggplot()+
  geom_line(data=ts_df, aes(x=time, y=ts, colour="1"), size=1.5)+
  geom_ribbon(data=ts_df, aes(x=year, y=ts, ymin=low_sd, ymax=high_sd), alpha=0.2)+
  geom_line(data=cod_ts_df, aes(x=time, y=ts, colour="2"), size=1.4, linetype="dashed")+
  geom_line(data=ampl_ts_df, aes(x=time, y=ts, colour="3"), size=1.4, linetype="dashed")+
  scale_color_manual(name="Both Species Model",values=c("#656565","#035d8a","#f4a582"),
                     labels=c("Estimated N","Atlantic cod","American plaice"))+
  ggtitle("Both species")+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text=element_text(size=16),legend.key = element_rect(fill = NA),legend.title=element_text(size=16, face="bold"),
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))

# create function to extract legend only from plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

#extract legend for first plot
legend <- get_legend(ampl_curve)

jpeg("functional_response_curves_w_error.jpeg", width=12, height=10, units = "in", res=400)
gridExtra::grid.arrange(arrangeGrob(ampl_curve+ theme(legend.position="none"),legend, nrow=1), 
                        arrangeGrob(cod_curve, both_curve, nrow=1), nrow=2,
                        bottom = textGrob("Relative Prey Density", vjust = 0, gp = gpar(fontface = "bold", cex = 1.5)),
                        left = textGrob("Probability of Consumption", rot = 90, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
dev.off()

#extract both legends for second plot
legend <- get_legend(ampl_ts)
legend2 <- get_legend(both_ts)

library(cowplot)
leg12<-plot_grid(legend, legend2, ncol=1)

jpeg("sand_lance_ts_zeroes_7_30.jpeg", width=12, height=10, units = "in", res=400)
gridExtra::grid.arrange(arrangeGrob(ampl_ts+ theme(legend.position="none"), leg12, nrow=1), 
                        arrangeGrob(cod_ts,both_ts+ theme(legend.position="none"), nrow=1), nrow=2,
                        bottom = textGrob("Year", vjust = 0, gp = gpar(fontface = "bold", cex = 1.5)),
                        left = textGrob("Sand Lance Abundance Index", rot = 90, vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)))
dev.off()


