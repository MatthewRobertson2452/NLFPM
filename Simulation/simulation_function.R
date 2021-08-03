

#Create simulation function

simulate_sampling<-function(
  nyears = 25,
  predator_sampling_type = "type1", #options are "random", "type1", "type2", and "type3"
  trawl_sampling = "binomial", #options are "ZIP" and "binomial"
  prey_trend = "RW", #options are "RW", steady" and "2states". 2states currently requires that nyears is divisible by 2
  RW_var_prey = 200, #variability in the normal distribution for the prey Random Walk
  pred_trend = "RW", #options are "RW", "steady" and "increasing" 
  RW_var_pred = 50, #variability in the normal distribution for the predator Random Walk
  prey_mean = 10000, #size of prey population. If two states this should be a vector of length 2
  pred_mean = 1000,  #size of the predation population. If increasing need a vector of length 2, first value is min second is max
  n_samples = 20, #number of trawl samples
  n_pred_samples = 20, #number of trawl samples of predators
  max_pred_tow = 20, #maximum number of predators whose stomachs are sampled in a tow
  min_pred_tow = 10, #minimum number of predators whose stomachs are sampled in a tow
  pred_obs_error = FALSE, #not functional yet
  attack = 10, #attack rate and handling rate for type 2 functional response if using
  handling= 1,
  intercept=c(0, 0),
  max_val=c(0.8,0.8),
  K=0.95,
  chi=0.3,
  beta=3)
{ 
  
  ## Prey population size
  #STEADY
  if(prey_trend=="steady"){
    prey_num<-rpois(n=nyears, prey_mean) #population size stationary around a mean
  }
  
  #2 STATES
  if(prey_trend=="2states"){
    prey_num<-c(rpois(n=nyears/2, prey_mean[1]),rpois(n=nyears/2, prey_mean[2])) #population stationary around two means at separate times
  }
  
  #RANDOM WALK
  if(prey_trend=="RW"){
    prey_num<-rep(NA, nyears)
    prey_num[1]<-rnorm(1, prey_mean, 50)
    for(i in 2:nyears){
      prey_num[i]<-prey_num[i-1]+rnorm(1, 0, RW_var_prey)
    }
    prey_num<-as.integer(prey_num)
  }

    ## Predator population size
  #STEADY
  if(pred_trend == "steady"){
    pred_num<-rpois(n=nyears, pred_mean) #population size stationary around a mean
  }
  
  #INCREASING
  if(pred_trend == "increasing"){
    pred_slope<-(pred_mean[2]-pred_mean[1])/nyears 
    years<-seq(from=0, to=24, by=1)
    pred_num<-rpois(nyears, pred_mean[1]+pred_slope*years)
  }
  
  #RANDOM WALK
  if(pred_trend == "RW"){
    pred_num<-rep(NA, nyears)
    pred_num[1]<-rnorm(1, pred_mean, 50)
    for(i in 2:nyears){
      pred_num[i]<-pred_num[i-1]+rnorm(1, 0, RW_var_pred)
    }
    pred_num<-as.integer(pred_num)
  }
  
  #Trawl sampling
  #can assume that samples are randomly drawn from the same mean as the population for now, may make more sense to scale it to a smaller value (catchability)
  #ZIP
  if(trawl_sampling=="ZIP"){
    std_trawl <- (prey_num - mean(prey_num))/sd(prey_num)
    trawl_prop <- exp(std_trawl)/(1+exp(std_trawl))
    
    #samples drawn from ZIP with equal mean through time but the probability of success is proportional to size of population
    #since I plan to treat trawl data as prsence/absence all values are then turned into binary
    bottom_trawl_df<-data.frame(year = rep(1, n_samples), n_pred = seq(from=1, to=n_samples, by=1), sl_presence = ifelse(rZIP(n_samples, 50, 1-trawl_prop[1])>0,1,0))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, n_samples), n_pred = seq(from=1, to=n_samples, by=1), sl_presence = ifelse(rZIP(n_samples, 50, 1-trawl_prop[i])>0,1,0))
      bottom_trawl_df<-rbind(bottom_trawl_df, year_one_df)
    }
  }
  
  #BINOMIAL
  if(trawl_sampling=="binomial"){
    std_trawl <- (prey_num - mean(prey_num))/sd(prey_num)
    logit_t <- exp(std_trawl)/(1+exp(std_trawl))
    
    trawl_prop<-1-exp(-logit_t)
    
    #samples drawn from binomial distribution where the probability of success is proportional to size of population
    bottom_trawl_df<-data.frame(year = rep(1, n_samples), n_pred = seq(from=1, to=n_samples, by=1), sl_presence = rbinom(n_samples,1, trawl_prop[1]))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, n_samples), n_pred = seq(from=1, to=n_samples, by=1), sl_presence = rbinom(n_samples,1, trawl_prop[i]))
      bottom_trawl_df<-rbind(bottom_trawl_df, year_one_df)
    }
  }
  
  ## Predation
  #random predation, does not vary with the number of prey
  if(predator_sampling_type=="random"){
    predation_df<-data.frame(year = rep(1, pred_num[1]), n_pred = seq(from=1, to=pred_num[1], by=1), sl_presence = rbinom(pred_num[1], 1, 0.5))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, pred_num[i]), n_pred = seq(from=1, to=pred_num[i], by=1), sl_presence = rbinom(pred_num[i], 1, 0.5))
      predation_df<-rbind(predation_df,year_one_df)
    }
  }
  
  #type 1 functional response, consumption probability depends on number of prey
  #minimum value of 0.1 and denominator multiplied by 2 so it doesn't scale fully bwn 0 and 1
  if(predator_sampling_type=="type1"){
    std_type1 <- (prey_num - mean(prey_num))/sd(prey_num)
    type_1_prop <- exp(std_type1)/(1+exp(std_type1))
    
    predation_df<-data.frame(year = rep(1, pred_num[1]), n_pred = seq(from=1, to=pred_num[1], by=1), sl_presence = rbinom(pred_num[1], 1, type_1_prop[1]))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, pred_num[i]), n_pred = seq(from=1, to=pred_num[i], by=1), sl_presence = rbinom(pred_num[i], 1, type_1_prop[i]))
      predation_df<-rbind(predation_df,year_one_df)
    }
  }
  
  #type 2 functional response, consumption probability depends on number of prey
  #minimum value of 0.1 and denominator multiplied by 2 so it doesn't scale fully bwn 0 and 1
  if(predator_sampling_type=="type2"){
    
    type_2<-((attack*prey_num)/(1+attack*handling*prey_num))
    type_2_sc <- (type_2 - mean(type_2))/sd(type_2)
    type_2_prop <- exp(type_2_sc)/(1+exp(type_2_sc))
    
    predation_df<-data.frame(year = rep(1, pred_num[1]), n_pred = seq(from=1, to=pred_num[1], by=1), sl_presence = rbinom(pred_num[1], 1, type_2_prop[1]))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, pred_num[i]), n_pred = seq(from=1, to=pred_num[i], by=1), sl_presence = rbinom(pred_num[i], 1, type_2_prop[i]))
      predation_df<-rbind(predation_df,year_one_df)
    }
  }
  
  if(predator_sampling_type=="type3"){
    
    std_trawl <- (prey_num - mean(prey_num))/sd(prey_num)
    logit_t <- exp(std_trawl)/(1+exp(std_trawl))
    
    trawl_prop<-1-exp(-logit_t)
    type_3_prop <- ((K*(trawl_prop^beta))/((chi^beta)+(trawl_prop^beta)))
    
    predation_df<-data.frame(year = rep(1, pred_num[1]), n_pred = seq(from=1, to=pred_num[1], by=1), sl_presence = rbinom(pred_num[1], 1, type_3_prop[1]))
    for(i in 2:nyears){
      year_one_df<-data.frame(year = rep(i, pred_num[i]), n_pred = seq(from=1, to=pred_num[i], by=1), sl_presence = rbinom(pred_num[i], 1, type_3_prop[i]))
      predation_df<-rbind(predation_df,year_one_df)
    }
  }
  
  ## Bottom trawl and stomach sampling of predators
  trawl_pred<-  round((pred_num - min(pred_num))/(max(pred_num)-min(pred_num)) * ((max_pred_tow-min_pred_tow)+min_pred_tow))
  
  npred_catch<-rpois(n_pred_samples, trawl_pred[1]) #how many predators were caught in the trawl survey
  npred_catch<-ifelse(npred_catch==0,1,npred_catch) #ensures that one predator is caught on each, could remove but easier for plotting now
  
  one_catch<-subset(predation_df, year==1) #initiate loop
  caught_fish<-sample(one_catch$n_pred, sum(npred_catch), replace=FALSE) #sample stomachs from the predators that were caught without replacement
  one_catch$pred_catch<-0
  one_catch$pred_catch[caught_fish]<-1 #whichever fish had their stomachs sampled give a value of 1
  pred_catch<-one_catch
  
  for(i in 2:nyears){
    npred_catch<-rpois(n_pred_samples, trawl_pred[i])
    npred_catch<-ifelse(npred_catch==0,1,npred_catch)
    one_catch<-subset(predation_df, year==i)
    caught_fish<-sample(one_catch$n_pred, sum(npred_catch), replace=FALSE)
    one_catch$pred_catch<-0
    one_catch$pred_catch[caught_fish]<-1
    pred_catch<-rbind(pred_catch, one_catch)
  }
  
  pred_catch_df<-subset(pred_catch, pred_catch>0) #create a dataframe for all predators that had their stomachs sampled
  pred_catch_df$obs_sl<-pred_catch_df$sl_presence
  if(pred_obs_error==TRUE){
    for(i in 1:length(pred_catch_df$year)){
      pred_catch_df$obs_sl[i]<-ifelse(pred_catch_df$obs_sl[i]==1,rbinom(1,1,1-obs_error),0)
    }
  }
  
  
  output_list<-list(
    bottom_trawl_df = bottom_trawl_df,
    pred_catch_df = pred_catch_df,
    prey_num = prey_num, 
    pred_num = pred_num,
    all_preds = pred_catch #all predators in the population
  )
  return(output_list)
}
