setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting_F.R")



#read data from JHU
JHU_F_data <- read.csv("JHU_data/time_series_covid19_deaths_global.csv")
JHU_F_data<-t(JHU_F_data)
JHU_C_data<-read.csv("JHU_data/time_series_covid19_confirmed_global.csv")
JHU_C_data<-t(JHU_C_data)
JHU_R_data<-read.csv("JHU_data/time_series_covid19_recovered_global.csv")
JHU_R_data<-t(JHU_R_data)

#this data goes up to May 4, for comparison with Flaxman et. al
C_data_May_4<-read.csv("JHU_data/time_series_covid19_confirmed_global - to May 4.csv")
C_data_May_4<-t(C_data_May_4)

F_data_May_4 <- read.csv("JHU_data/time_series_covid19_deaths_global - to May 4.csv")
F_data_May_4<-t(F_data_May_4)

#generates a vector of region labels for the estimates (JHU data)
regions=function(data){
  
  regions=c()
  for(x in seq(1, ncol(data))){
    #combines label from first and second row (region and country)
    regions=append(regions, paste(data[1, x], data[2, x]))
  }
  return(regions)
}

#fits r_tilde, p, alpha, and K_tilde for a given region
fit_param_F=function(region, data){

  #format case data for a given region
  cases=as.integer(data[5:nrow(data), region])
  cases=cases[!is.na(cases)]
  #discards data for days before initial outbreak
  start=min(which(cases>0, arr.ind=TRUE))
  cases=c(cases[start:length(cases)])
  
  #starting guess for parameters
  par=c(r_tilde=2, p=1, alpha=1, K_tilde=cases[length(cases)])
  
  #performs the fit
  fit=optim(par=par, fn=ssq_F, cases=cases, control=list(parscale=c(1,1,1,10^floor(log10(cases[length(cases)])))))
  parest=fit$par
  
  return(parest)
}

fit_multiple_F=function(data){
  
  #removes columns for countries with no cases
  few_cases=c()
  for(x in seq(1, ncol(data))){
    cases=as.integer(data[5:nrow(data), x])
    cases=cases[!is.na(cases)]
    if(cases[length(cases)]<=10){
      few_cases=append(few_cases, x)
    }
  }
  data=subset(data, select=-c(few_cases))
  
  #generates labels
  regions=regions(data)
  
  #performs the fit for the first element in the dataframe (assumption: the dataframe has at least 2 columns)
  paramdf=data.frame(fit_param_F(1, data))
  colnames(paramdf)=c(regions[1])
  
  #performs the fit for the rest of the regions, adds results to the dataframe
  #if all the data were working, the loop would be for seq(2, ncol(data))
  for(x in seq(2, ncol(data))){
    print(x)
    p=data.frame(fit_param_F(x, data))
    colnames(p)=c(regions[x])
    paramdf=cbind(paramdf, p)
    print(p)
  }
  return(paramdf)
}

fit_tau_mu_CFR=function(region, C_data=JHU_C_data, F_data=JHU_F_data){
  
  #estimates paramters for fatality data (discarding days before outbreak)
  F_parest=fit_param_F(region, F_data)
  
  #F_parest=c(F_parest[1], F_parest[2], F_parest[3], F_parest[4])
  
  #format case data for a given region
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  #format fatality data for a given region
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]

  #use this if just fitting tau
  #par=c(tau=20)
  # fit=optim(par=par, fn=ssq_C_F, cases_C=cases_C, cases_F=cases_F, F_parest=F_parest, method="Brent", lower=-30, upper=30)
  
  #for fitting tau and mu_CFR
  par=c(tau=20, mu_CFR=0.05)
  fit=optim(par=par, fn=ssq_C_F, cases_C=cases_C, cases_F=cases_F, F_parest=F_parest, control=list(parscale=c(1,.5)))
  
  parest=fit$par
  
  return(parest)
}

#plots the unfitted case data scaled down to visualize the curves in comparison to one another
plot_cases_scaled=function(region, C_data=JHU_C_data, F_data=JHU_F_data, factor=0.05){
  
  regions=regions(C_data)
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  
  times=c(1:length(cases_C))
  
  C_df=data.frame(times, cases_C)
  F_df=data.frame(times, cases_F)
  
  plot=ggplot(data=C_df, aes(x=times, y=cases_C*factor, color="green"))+geom_line()+geom_line(data=F_df, aes(x=times, y=cases_F, color="red"))+theme(legend.position="none")+labs(title=regions[region])
  print(plot)
}

#plots the shifted and scaled fitted curves to see how good the tau and mu_CFR fit is
#par is returned from fit_tau_mu_CFR
plot_shifted_scaled_cases=function(region, par=NULL, C_data=JHU_C_data, F_data=JHU_F_data){
  if(is.null(par)==TRUE){
    par=fit_tau_mu_CFR(region, C_data, F_data)
    }

  
  regions=regions(C_data)
  
  tau=par[1]
  factor=par[2]
  
  #fits parameters to generate values for F
  F_parest=fit_param_F(region, F_data)
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  
  #shift times for cumulative cases by tau
  times_C=c(1:length(cases_C))+tau
  #C_df=data.frame(times_C, cases_C)
  
  start=min(which(cases_F>0, arr.ind=TRUE))
  init=cases_F[start]
  
  #generates values for F from fit
  F_df=generate_F(times_C, init, start, F_parest)
  
  #generates values for C from fit
  par_C=c(r=2, p=1, alpha=1, K=5000)
  C_df=generate_C(par_C, cases_C, times_C)
  
  plot=ggplot(data=C_df, aes(x=times, y=cases_C*factor, color="red"))+geom_line()+geom_line(data=F_df, aes(x=times, y=y, color="green"))+geom_line()+theme(legend.position="none")+labs(title=regions[region])
  print(plot)
}

#for this function, factor is the estimate for mu_CFR returned from fit_tau_mu_CFR, and mu_CFR is a sequence of guesses for the true mu_CFR
#function plots the proportion of cases reported vs. actual mu_CFR
plot_underreporting_vs_mu_CFR=function(factor, mu_CFR, region, C_data=JHU_C_data, F_data=JHU_F_data){
  
  regions=regions(C_data)
  ratios=c()
  
  for(x in mu_CFR){
    ratio=x/factor
    ratios=append(ratios, ratio)
  }
  
  df=data.frame(mu_CFR, ratios)
  
  #print(ratios[1])
  #print(ratios[length(ratios)])
  
  plot=ggplot(data=df, aes(x=mu_CFR, y=ratios))+geom_line()+ labs(title=regions[region])
  print(plot)
}

#returns the "underreporting ratio" vs time
phi_vs_time=function(region, C_data=JHU_C_data, F_data=JHU_F_data, mu_CFR=0.01){

  parest=fit_tau_mu_CFR(region, C_data, F_data)
  
  tau=parest[1]
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  
  start=min(which(cases_F>0, arr.ind=TRUE))
  init=cases_F[start]
  F_parest=fit_param_F(region, F_data)
  
  #shift times for cumulative cases by tau
  times=c(1:length(cases_C))+tau
 
  #generates C and F curves from the fits
  par_C=c(r=2, p=1, alpha=1, K=5000)
  C_df=generate_C(par_C, cases_C, times)
  
  F_df=generate_F(times, init, start, F_parest)
  
  #starts calculating the ratio at the point where both C and F become nonzero (should be very close together, since they have been shifted)
  #use the fit instead of the data?
  start=max(min(which(cases_C>0, arr.ind=TRUE)), min(which(cases_F>0, arr.ind=TRUE)))
  
  #removes the zeroes before "start"
  F_df=F_df[-(1:start-1),]
  C_df=C_df[-(1:start-1),]
  times=times[-(1:start-1)]
  
  ratios=((1/mu_CFR)*F_df$y)/C_df$cases_C

  df=data.frame(times, ratios)
  return(df)
}

#plots the ratios calculated in ratios_vs_time()
plot_phi_vs_time=function(region, C_data=JHU_C_data, F_data=JHU_F_data, mu_CFR=.01){
  regions=regions(C_data)
  
  df=phi_vs_time(region, C_data, F_data, mu_CFR=mu_CFR)
  plot=ggplot(data=df, aes(x=times, y=ratios))+geom_line()+labs(title=regions[region])+theme(legend.position="none")
  print(plot)
}

#returns an estimate for total number of people infected in a region to date, accounting for underreporting
total_infected=function(region, C_data=JHU_C_data, F_data=JHU_F_data, mu_CFR=0.01){
  regions=regions(C_data)
  
  df=phi_vs_time(region, C_data, F_data, mu_CFR=mu_CFR)
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  
  #starts calculating the ratio at the point where both C and F become nonzero (should be very close together, since they have been shifted)
  start=max(min(which(cases_C>0, arr.ind=TRUE)), min(which(cases_F>0, arr.ind=TRUE)))
  pre_start_cases=cases_C[start-1]
  cases_C=c(cases_C[start:length(cases_C)])
  
  
  new_cases=c(cases_C[1])
  for(x in c(2:length(cases_C))){
    new_case=cases_C[x]-cases_C[x-1]
    new_cases=append(new_cases, new_case)
  }
  
  true_new_cases=new_cases*df$ratios

  total_cases=sum(true_new_cases)+pre_start_cases
  
  print(regions[region])

  return(total_cases)
}

total_infected_multiple_regions=function(C_data=JHU_C_data, F_data=JHU_F_data, mu_CFR=0.01){
  #removes columns for countries with no cases
  few_cases=c()
  for(x in seq(1, ncol(C_data))){
    cases=as.integer(C_data[5:nrow(C_data), x])
    cases=cases[!is.na(cases)]
    if(cases[length(cases)]<=10){
      few_cases=append(few_cases, x)
    }
  }
  C_data=subset(C_data, select=-c(few_cases))
  F_data=subset(F_data, select=-c(few_cases))
  
  #generates labels
  regions=regions(C_data)
  
  totals=c()
  
  for(x in seq(1, ncol(C_data))){
    print(x)
    total=total_infected(region=x, C_data=C_data, F_data=F_data, mu_CFR=mu_CFR)
    totals=append(totals, total)
    print(total)
  }
  
  return(totals)
}
  
