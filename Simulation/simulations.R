
library(TMB)

#prep models so theyll run
dyn.unload("LFPM")
compile("LFPM.cpp")
dyn.load("LFPM")

dyn.unload("NLFPM")
compile("NLFPM.cpp")
dyn.load("NLFPM")

dyn.unload("trawl_only_model")
compile("trawl_only_model.cpp")
dyn.load("trawl_only_model")

###############################################
#######################################
############TYPE 2 SIMULATIONS
#######################################
#############################################

#define number of simulations
nsim<-1000

#create storage objects for everything we're interested in
mse_diff<-matrix(NA, nrow=nsim, ncol=3)
aic_diff<-matrix(NA, nrow=nsim, ncol=2)
true_mse<-matrix(NA, nrow=nsim, ncol=3)
iye_bias<-matrix(NA, nrow=nsim, ncol=3)
hess_mat<-matrix(NA, nrow=nsim, ncol=1)
grad_mat<-matrix(NA, nrow=nsim, ncol=1)
sel_k_mat<-matrix(NA, nrow=nsim, ncol=1)
sel_chi_mat<-matrix(NA, nrow=nsim, ncol=1)
beta_mat<-matrix(NA, nrow=nsim, ncol=1)
sd_diff<-matrix(NA, nrow=50, ncol=nsim)

#run the simulations
for(k in 1:nsim){
  
  print(rep(k,50)) #print to keep track of what run we are run
  
  #run simulation function
  sample_dat<-simulate_sampling(nyears=25, max_pred_tow=20, min_pred_tow=10,  n_samples = 200, #number of trawl samples
                                n_pred_samples = 100, pred_mean=10000, RW_var_pred=200,
                                predator_sampling_type = "type3", trawl_sampling="binomial", K=1,
                                chi=0.3, beta=1)
  
  #mean trawl proportion
  trawl_agg<-aggregate(sample_dat$bottom_trawl_df$sl_presence, by=list(sample_dat$bottom_trawl_df$year), FUN=mean)
  #mean predator proportion
  pred_agg<-aggregate(sample_dat$pred_catch_df$obs_sl, by=list(sample_dat$pred_catch_df$year), FUN=mean)
  
  #prepare data to be used in the models
  nyrs<-length(unique(sample_dat$bottom_trawl_df$year))
  pa<-c(sample_dat$bottom_trawl_df$sl_presence)
  n<-length(pa)
  ndex<-1
  iyear <- c(sample_dat$bottom_trawl_df$year-1)
  idex <- c(rep(0, length(sample_dat$bottom_trawl_df$year)))
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=1,
    iyear=iyear,
    idex=idex,
    k=1,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    logrw_var = log(2)
  )
  
  ###RUNNING A TYPE 1 MODEL
  obj <- MakeADFun(tmb.data,parameters,DLL="trawl_only_model",inner.control=list(maxit=2000,trace=F), random="iye", silent=TRUE)
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  nosl_rep = obj$report()
  
  nosl_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  sdrep_nosl<-sdreport(obj)
  
  nyrs<-length(unique(sample_dat$bottom_trawl_df$year))
  pa<-c(sample_dat$bottom_trawl_df$sl_presence, sample_dat$pred_catch_df$obs_sl)
  n<-length(pa)
  ndex<-2
  iyear <- c(sample_dat$bottom_trawl_df$year-1, sample_dat$pred_catch_df$year-1)
  idex <- c(rep(0, length(sample_dat$bottom_trawl_df$year)), rep(1, length(sample_dat$pred_catch_df$year)))
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=2,
    iyear=iyear,
    idex=idex,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    lk=0
  )
  
  obj <- MakeADFun(tmb.data,parameters,DLL="LFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  t1_rep = obj$report()
  
  t1_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=2,
    iyear=iyear,
    idex=idex,
    k=1,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    lbeta = log(2),
    lchi = log(0.5),
    logrw_var=log(2)
  )
  
  
  obj <- MakeADFun(tmb.data,parameters,DLL="NLFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
  
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  t2_rep = obj$report()
  
  sdrep<-sdreport(obj)
  
  sd_diff[,k]<-sdrep$sd-sdrep_nosl$sd
  
  grad_mat[k,1]<-max(obj$gr())
  
  t2_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  sel_chi_mat[k,1]<-t2_rep$chi
  
  hess_mat[k,1]<-sdrep$pdHess
  
  beta_mat[k,1]<-t2_rep$beta
  
  sc_nosl<-nosl_rep$iye
  sc_t1<-t1_rep$iye
  sc_tau<-t2_rep$iye
  
  std_trawl <- (sample_dat$prey_num - mean(sample_dat$prey_num))/sd(sample_dat$prey_num)
  sc_prey <- exp(std_trawl)/(1+exp(std_trawl))

  
  mse_diff[k,1]<-mean((sc_prey-sc_tau)^2)-mean((sc_prey-sc_t1)^2)
  mse_diff[k,2]<-mean((sc_prey-sc_tau)^2)-mean((sc_prey-sc_nosl)^2)
  mse_diff[k,3]<-mean((sc_prey-sc_t1)^2)-mean((sc_prey-sc_nosl)^2)
  
  aic_diff[k,1]<-t2_aic-t1_aic
  
  true_mse[k,1]<-mean((sc_prey-sc_tau)^2)
  true_mse[k,2]<-mean((sc_prey-sc_t1)^2)
  true_mse[k,3]<-mean((sc_prey-sc_nosl)^2)
  
  iye_bias[k,1]<-mean(sc_tau-sc_prey)
  iye_bias[k,2]<-mean(sc_t1-sc_prey)
  iye_bias[k,3]<-mean(sc_nosl-sc_prey)
}


#create list to save all output objects
sim_run_type2<-list(
  mse_diff=mse_diff,
  aic_diff=aic_diff,
  true_mse=true_mse,
  iye_bias=iye_bias,
  hess_mat=hess_mat,
  grad_mat=grad_mat,
  chi_mat=sel_chi_mat,
  beta_mat=beta_mat,
  sd_diff=sd_diff
)  

save(sim_run_type2, file="sim_run_type2_w_nosl_mod_7_30.RData")

###############################################
#######################################
############TYPE 2 SIMULATIONS
#######################################
#############################################

#define number of simulations
nsim<-1000

#create storage objects for everything we're interested in
mse_diff<-matrix(NA, nrow=nsim, ncol=3)
aic_diff<-matrix(NA, nrow=nsim, ncol=2)
true_mse<-matrix(NA, nrow=nsim, ncol=3)
iye_bias<-matrix(NA, nrow=nsim, ncol=3)
hess_mat<-matrix(NA, nrow=nsim, ncol=1)
grad_mat<-matrix(NA, nrow=nsim, ncol=1)
sel_k_mat<-matrix(NA, nrow=nsim, ncol=1)
sel_chi_mat<-matrix(NA, nrow=nsim, ncol=1)
beta_mat<-matrix(NA, nrow=nsim, ncol=1)
sd_diff<-matrix(NA, nrow=50, ncol=nsim)

#run the simulations
for(k in 1:nsim){
  
  print(rep(k,50)) #print to keep track of what run we are run
  
  #run simulation function
  sample_dat<-simulate_sampling(nyears=25, max_pred_tow=20, min_pred_tow=10,  n_samples = 200, #number of trawl samples
                                n_pred_samples = 100, pred_mean=10000, RW_var_pred=200,
                                predator_sampling_type = "type3", trawl_sampling="binomial", K=1,
                                chi=0.3, beta=3)
  
  #mean trawl proportion
  trawl_agg<-aggregate(sample_dat$bottom_trawl_df$sl_presence, by=list(sample_dat$bottom_trawl_df$year), FUN=mean)
  #mean predator proportion
  pred_agg<-aggregate(sample_dat$pred_catch_df$obs_sl, by=list(sample_dat$pred_catch_df$year), FUN=mean)
  
  #prepare data to be used in the models
  nyrs<-length(unique(sample_dat$bottom_trawl_df$year))
  pa<-c(sample_dat$bottom_trawl_df$sl_presence)
  n<-length(pa)
  ndex<-1
  iyear <- c(sample_dat$bottom_trawl_df$year-1)
  idex <- c(rep(0, length(sample_dat$bottom_trawl_df$year)))
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=1,
    iyear=iyear,
    idex=idex,
    k=1,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    logrw_var = log(2)
  )
  
  ###RUNNING A TYPE 1 MODEL
  obj <- MakeADFun(tmb.data,parameters,DLL="trawl_only_model",inner.control=list(maxit=2000,trace=F), random="iye", silent=TRUE)
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  nosl_rep = obj$report()
  
  nosl_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  sdrep_nosl<-sdreport(obj)
  
  nyrs<-length(unique(sample_dat$bottom_trawl_df$year))
  pa<-c(sample_dat$bottom_trawl_df$sl_presence, sample_dat$pred_catch_df$obs_sl)
  n<-length(pa)
  ndex<-2
  iyear <- c(sample_dat$bottom_trawl_df$year-1, sample_dat$pred_catch_df$year-1)
  idex <- c(rep(0, length(sample_dat$bottom_trawl_df$year)), rep(1, length(sample_dat$pred_catch_df$year)))
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=2,
    iyear=iyear,
    idex=idex,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    lk=0
  )
  
  obj <- MakeADFun(tmb.data,parameters,DLL="LFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  t1_rep = obj$report()
  
  t1_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  
  tmb.data<-list(
    n=n,
    nyrs=nyrs,
    ndex=2,
    iyear=iyear,
    idex=idex,
    k=1,
    pa=pa)
  
  parameters<-list(
    iye = seq(from=log(10), to=log(20), length.out=tmb.data$nyrs),
    lbeta = log(2),
    lchi = log(0.5),
    logrw_var=log(2)
  )
  
  
  obj <- MakeADFun(tmb.data,parameters,DLL="NLFPM",inner.control=list(maxit=2000,trace=F), random=c("iye"), silent=TRUE)
  
  opt<-nlminb(obj$par,obj$fn,obj$gr,control = list(trace=10,eval.max=2000,iter.max=1000),silent=TRUE)
  t2_rep = obj$report()
  
  sdrep<-sdreport(obj)
  
  sd_diff[,k]<-sdrep$sd-sdrep_nosl$sd
  
  grad_mat[k,1]<-max(obj$gr())
  
  t2_aic<-2*opt$objective+ 2*length(opt$par)#aic
  
  sel_chi_mat[k,1]<-t2_rep$chi
  
  hess_mat[k,1]<-sdrep$pdHess
  
  beta_mat[k,1]<-t2_rep$beta
  
  sc_nosl<-nosl_rep$iye
  sc_t1<-t1_rep$iye
  sc_tau<-t2_rep$iye
  
  std_trawl <- (sample_dat$prey_num - mean(sample_dat$prey_num))/sd(sample_dat$prey_num)
  sc_prey <- exp(std_trawl)/(1+exp(std_trawl))
  
  
  mse_diff[k,1]<-mean((sc_prey-sc_tau)^2)-mean((sc_prey-sc_t1)^2)
  mse_diff[k,2]<-mean((sc_prey-sc_tau)^2)-mean((sc_prey-sc_nosl)^2)
  mse_diff[k,3]<-mean((sc_prey-sc_t1)^2)-mean((sc_prey-sc_nosl)^2)
  
  aic_diff[k,1]<-t2_aic-t1_aic
  
  true_mse[k,1]<-mean((sc_prey-sc_tau)^2)
  true_mse[k,2]<-mean((sc_prey-sc_t1)^2)
  true_mse[k,3]<-mean((sc_prey-sc_nosl)^2)
  
  iye_bias[k,1]<-mean(sc_tau-sc_prey)
  iye_bias[k,2]<-mean(sc_t1-sc_prey)
  iye_bias[k,3]<-mean(sc_nosl-sc_prey)
}


#create list to save all output objects
sim_run_type3<-list(
  mse_diff=mse_diff,
  aic_diff=aic_diff,
  true_mse=true_mse,
  iye_bias=iye_bias,
  hess_mat=hess_mat,
  grad_mat=grad_mat,
  chi_mat=sel_chi_mat,
  beta_mat=beta_mat,
  sd_diff=sd_diff
)  

save(sim_run_type3, file="sim_run_type3_w_nosl_mod_7_30.RData")


load("sim_run_type3_w_nosl_mod_7_30.RData")
load("sim_run_type2_w_nosl_mod_7_30.RData")

nsim<-1000

type2<-sim_run_type2
type3<-sim_run_type3

library(ggplot2)
library(reshape2)

t2_mse_df<-data.frame(mse_diff=c(type2$true_mse[,3],type2$true_mse[,2],type2$true_mse[,1]), test=c(rep("Trawl", nsim),rep("LFPM", nsim), rep("NLFPM", nsim)))
t3_mse_df<-data.frame(mse_diff=c(type3$true_mse[,3],type3$true_mse[,2],type3$true_mse[,1]), test=c(rep("Trawl", nsim),rep("LFPM", nsim), rep("NLFPM", nsim)))

jpeg("mse_type2_3_wtrawl_7_30.jpeg", units="in", width=10, height=5, res=400)
p1<-ggplot(t2_mse_df, aes(x=test, y=mse_diff, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#035d8a","#f4a582","#c2a5cf"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Estimation Model")+ylab("Mean Squared Error")+ggtitle("Type II Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))

p2<-ggplot(t3_mse_df, aes(x=test, y=mse_diff, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+scale_fill_manual(values=c("#035d8a","#f4a582","#c2a5cf"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Estimation Model")+ylab("Mean Squared Error")+ggtitle("Type III Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=16,face="bold"))

gridExtra::grid.arrange(p1,p2, nrow=1)
dev.off()

aic_diff_df<-data.frame(aic_diff=c(type2$aic_diff[,1],type3$aic_diff[,1]), test=c(rep("Type II", nsim), rep("Type III", nsim)))

jpeg("aic_diff_type2_3_7_30.jpeg", units="in", width=5, height=5, res=400)
ggplot(aic_diff_df, aes(x=test, y=aic_diff, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  ylim(-1700, 10)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de","#f4a582"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Operating Model")+ylab("AIC Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"))

dev.off()


t2_bias_df<-data.frame(bias=c(type2$iye_bias[,1],type2$iye_bias[,2],type2$iye_bias[,3]), test=c(rep("NLFPM", nsim), rep("LFPM", nsim), rep("Trawl", nsim)))

t3_bias_df<-data.frame(bias=c(type3$iye_bias[,1],type3$iye_bias[,2], type3$iye_bias[,3]), test=c(rep("NLFPM", nsim), rep("LFPM", nsim), rep("Trawl", nsim)))


jpeg("bias_type2_3_wtrawl_7_30.jpeg", units="in", width=10, height=5, res=400)
p1<-ggplot(t2_bias_df, aes(x=test, y=bias, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#035d8a","#f4a582","#c2a5cf"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Estimation Model")+ylab("Bias")+ggtitle("Type II Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"), plot.title=element_text(size=16,face="bold"))

p2<-ggplot(t3_bias_df, aes(x=test, y=bias, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#035d8a","#f4a582","#c2a5cf"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Estimation Model")+ylab("Bias")+ggtitle("Type III Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"), plot.title=element_text(size=16,face="bold"))
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()





jpeg("parm_estimates_thicklines_7_30.jpeg", units="in", width=12, height=8, res=400)
chi_mat<-data.frame(chi=type2$chi_mat, parm=rep("chi",nsim))
p1<-ggplot(chi_mat, aes(x=parm, y=chi, fill=parm))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de"))+
  geom_hline(yintercept=0.3, linetype="dashed", colour="#ca0020", size=1)+
  xlab("Chi")+ylab("")+ylim(0.23,0.36)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0.5,1), "cm"))

beta_mat<-data.frame(beta=type2$beta_mat, parm=rep("chi",nsim))
p2<-ggplot(beta_mat, aes(x=parm, y=beta, fill=parm))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de"))+
  geom_hline(yintercept=1, linetype="dashed", colour="#ca0020", size=1)+
  xlab("Beta")+ylab("")+ylim(0.7,1.25)+ggtitle("Type II Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold",hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0.5,1), "cm"))

new_x<-seq(from=0, to=1, by=0.01)
K<-1
chi<-0.3
beta<-1
type3_true<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
K<-median(rep(1,nsim), na.rm=TRUE)
chi<-median(type2$chi_mat, na.rm=TRUE)
beta<-median(type2$beta_mat, na.rm=TRUE)
type3_model<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
curves_df<-data.frame(y=c(type3_true, type3_model), x=c(new_x,new_x), group=c(rep("True", length(new_x)),rep("Model", length(new_x))))

new_x_mat<-matrix(new_x, nrow=length(new_x), ncol=1000)
new_x_mat<-t(new_x_mat)

new_beta_mat<-matrix(type2$beta_mat, nrow=1000, ncol=length(new_x))
new_chi_mat<-matrix(type2$chi_mat, nrow=1000, ncol=length(new_x))

new_y_mat<-(1*new_x_mat^new_beta_mat)/((new_chi_mat^new_beta_mat)+(new_x_mat^new_beta_mat))

rownames(new_y_mat) = paste("trial", seq(1000), sep="")
colnames(new_y_mat) = paste("density", new_x, sep="")

dat = as.data.frame(new_y_mat)
dat$trial = rownames(dat)
mdat = melt(dat, id.vars="trial")
mdat$density = as.numeric(gsub("density", "", mdat$variable))

p3<-ggplot()+
  geom_line(data=mdat, aes(x=density, y=value, group=trial),size=1, alpha=0.03, col="grey")+
  geom_line(data=curves_df, aes(x=x, y=y, group=group, colour=group,linetype=group),size=1)+ 
  scale_color_manual(name="",labels=c("Median Estimated","True Dynamics"),values=c("#005cff", "#ca0020"))+
  scale_linetype_manual(name="",labels=c("Median Estimated","True Dynamics"),values=c("solid", "dashed"))+
  ylim(0,1)+
  xlab("Prey Density")+ylab("Probability of Consumption")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.8, 0.5),
        legend.text=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = NA),
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold"),
        plot.margin = unit(c(0,1,0,1), "cm"))



###############TYPE 3

chi_mat<-data.frame(chi=type3$chi_mat, parm=rep("chi",nsim))
p4<-ggplot(chi_mat, aes(x=parm, y=chi, fill=parm))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de"))+
  geom_hline(yintercept=0.3, linetype="dashed", colour="#ca0020", size=1)+
  xlab("Chi")+ylab("")+ylim(0.23,0.36)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0.5,1), "cm"))

beta_mat<-data.frame(beta=type3$beta_mat, parm=rep("chi",nsim))
p5<-ggplot(beta_mat, aes(x=parm, y=beta, fill=parm))+
  geom_jitter(color="darkgrey", size=0.5)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de"))+
  geom_hline(yintercept=3, linetype="dashed", colour="#ca0020", size=1)+
  xlab("Beta")+ylab("")+ylim(2.25,3.75)+ggtitle("Type III Operating Model")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold",hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,1,0.5,1), "cm"))

new_x<-seq(from=0, to=1, by=0.01)
K<-1
chi<-0.3
beta<-3
type3_true<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
K<-median(rep(1,nsim), na.rm=TRUE)
chi<-median(type3$chi_mat, na.rm=TRUE)
beta<-median(type3$beta_mat, na.rm=TRUE)
type3_model<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
curves_df<-data.frame(y=c(type3_true, type3_model), x=c(new_x,new_x), group=c(rep("True", length(new_x)),rep("Model", length(new_x))))


new_x_mat<-matrix(new_x, nrow=length(new_x), ncol=1000)
new_x_mat<-t(new_x_mat)

new_beta_mat<-matrix(type3$beta_mat, nrow=1000, ncol=length(new_x))
new_chi_mat<-matrix(type3$chi_mat, nrow=1000, ncol=length(new_x))

new_y_mat<-(1*new_x_mat^new_beta_mat)/((new_chi_mat^new_beta_mat)+(new_x_mat^new_beta_mat))

rownames(new_y_mat) = paste("trial", seq(1000), sep="")
colnames(new_y_mat) = paste("density", new_x, sep="")

dat = as.data.frame(new_y_mat)
dat$trial = rownames(dat)
mdat = melt(dat, id.vars="trial")
mdat$density = as.numeric(gsub("density", "", mdat$variable))

p6<-ggplot()+
  geom_line(data=mdat, aes(x=density, y=value, group=trial),size=1, alpha=0.03, col="grey")+
  geom_line(data=curves_df, aes(x=x, y=y, group=group, colour=group,linetype=group), size=1)+
  scale_color_manual(name="",labels=c("Median Estimated","True Dynamics"),values=c("#005cff", "#ca0020"))+
  scale_linetype_manual(name="",labels=c("Median Estimated","True Dynamics"),values=c("solid", "dashed"))+
  ylim(0,1)+
  xlab("Prey Density")+ylab("Probability of Consumption")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",
        legend.text=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = NA),
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold"),
        plot.margin = unit(c(1,1,0,1), "cm"))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

legend <- get_legend(p6)

gridExtra::grid.arrange( blankPlot, blankPlot, legend,
                         p1, p2, p3+theme(legend.position="none"),
                         p4, p5, p6+theme(legend.position="none"), nrow = 3, ncol=3,
                         heights=c(0.2,2.5,2.5))
dev.off()


sd_diff_df<-data.frame(sd_diff=c(c(type2$sd_diff),(type3$sd_diff)), test=c(rep("Type II", nsim*50), rep("Type III", nsim*50)))

jpeg("sd_diff_7_30.jpeg", units="in", width=5, height=5, res=400)
ggplot(sd_diff_df, aes(x=test, y=sd_diff, fill=test))+
  geom_jitter(color="darkgrey", size=0.5)+
  #ylim(-1700, 10)+
  geom_boxplot(alpha=0.5)+  scale_fill_manual(values=c("#92c5de","#f4a582"))+
  geom_hline(yintercept=0, linetype="dashed", colour="#ca0020", size=1.3)+
  xlab("Operating Model")+ylab("Standard Deviation Difference")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14,face="bold"))

dev.off()



new_x<-seq(from=0.01, to=1, by=0.01)
K<-1
chi<-0.3
beta<-3
type3_curve<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
K<-1
chi<-0.3
beta<-1
type2_curve<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
curves_df<-data.frame(y=c(type2_curve, type3_curve), x=c(new_x,new_x), group=c(rep("Type II", length(new_x)),rep("Type III", length(new_x))))

jpeg("estimated_shapes.jpeg", units="in", width=7, height=5, res=400)
ggplot(curves_df, aes(x=x, y=y, group=group, colour=group))+
  geom_line(size=1.4)+ scale_color_manual(name="",labels=c("Median Estimated","True Dynamics"),values=c("#92c5de","#f4a582"))+
  #ylim(0,1)+xlim(0,1)+
  xlab("Prey Density")+ylab("Probability of Consumption")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text=element_text(size=14),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text= element_text(size=14), axis.title=element_text(size=14), plot.title=element_text(size=16,face="bold"),
        plot.margin = unit(c(1,2,0,1), "cm"))+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip = 'off') + 
  annotate(geom="text",x=1.08, y=0.775, label="Type II", size=6,col="#92c5de")+
  annotate(geom="text",x=1.083, y=0.975, label="Type III", size=6,col="#f4a582")
dev.off()
