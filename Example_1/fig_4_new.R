# Plotting a scatter of posterior lines to represent the uncertainty.
#ploting the uncertinity of the posterior
# Plotting the data
par(mfrow = c(1,1))

plot(grnd.truth.I_, type="l", xlab = "Time(weeks)",
     ylab = "Proportion infected",ylim=c(0,0.35),col="white")

saved_estimations_mcmc=matrix(data=NA,ncol=length(grnd.truth.I_),nrow=length(MCMC.run[['gamma']])/step)

#MCMC
i=1
g=0
step=7 
while(i < length(MCMC.run[['gamma']])) {
  gamma_pre=(MCMC.run[['gamma']][i])
  beta_pre=(MCMC.run[['beta']][i])
  parms= c(beta= beta_pre,gamma=gamma_pre)
  Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
  Data= Data[,3]
  Data= Data[min_time:max_time]
  
  for(h in 1:1000){
    Data.noisy <- Data +
      rnorm(length(Data),0,sqrt(MCMC.run[['sigma2']][i]))
    if (sum(Data.noisy<0)==0){
      break
    }
  }
  
  saved_estimations_mcmc[g+1,]=Data.noisy
  g=g+1
  i=i+step
  print(i)
  #lines(true_time,Data.noisy,col="pink")
}

v1 <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
mynewdata_mcmc<-sapply(as.data.frame(saved_estimations_mcmc),function(i) quantile(i, v1))
mynewdata0_mcmc<-mynewdata_mcmc[1,]
mynewdata100_mcmc<-mynewdata_mcmc[7,]

polygon(c(1:length(grnd.truth.I_),length(grnd.truth.I_):1),c(mynewdata0_mcmc,rev(mynewdata100_mcmc)),col=rgb(1, 0, 0,0.3),border=NA)

lines(mynewdata0_mcmc,col="red")
lines(mynewdata100_mcmc,col="red")

#__________________________________________________________________
#ABC_SMC
i=1
g=0
step=5 
saved_estimations_smcabc_curr=matrix(data=NA,ncol=length(grnd.truth.I_),nrow=length(saved_beta_curr)/step)

while(i <= length(saved_beta_curr)) {
  gamma_pre=(saved_gamma_curr[i])
  beta_pre=(saved_beta_curr[i])
  parms= c(beta= beta_pre,gamma=gamma_pre)
  Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
  Data= Data[,3]
  Data= Data[min_time:max_time]
  
  for(h in 1:1000){
    Data.noisy <- Data +rnorm(length(Data),0,sqrt(saved_sigma2_curr[i]))
    if (sum(Data.noisy<0)==0){
      break
    }
  }
  saved_estimations_smcabc_curr[g+1,]=Data.noisy
  lines(Data.noisy)
  g=g+1
  
  i=i+step
  
  print(i)
  #lines(Data,col="red",lwd=2)
}

mynewdata_smcabc_curr<-sapply(as.data.frame(saved_estimations_smcabc_curr),function(i) quantile(i, v1))
mynewdata0_smcabc_curr<-mynewdata_smcabc_curr[1,]
mynewdata100_smcabc_curr<-mynewdata_smcabc_curr[7,]

polygon(c(1:length(grnd.truth.I_),length(grnd.truth.I_):1),c(mynewdata0_smcabc_curr,rev(mynewdata100_smcabc_curr)),col="green",border=NA)

lines(mynewdata0_smcabc_curr,col="green",lty=2)
lines(mynewdata100_smcabc_curr,col="green",lty=2)


#__________________________________________________________________
#ABC_SMC_new
i=1
g=0
step=1
saved_estimations_smcabc_modi=matrix(data=NA,ncol=length(grnd.truth.I_),nrow=length(saved_beta_modi)/step)

while(i <= length(saved_beta_modi)) {
  gamma_pre=(saved_gamma_modi[i])
  beta_pre=(saved_beta_modi[i])
  parms= c(beta= beta_pre,gamma=gamma_pre)
  Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
  Data= Data[,3]
  Data= Data[min_time:max_time]
  
  for(h in 1:1000){
    Data.noisy <- Data +rnorm(length(Data),0,sqrt(saved_sigma_modi[i]))
    if (sum(Data.noisy<0)==0){
      break
    }
  }
  saved_estimations_smcabc_modi[g+1,]=Data.noisy
  lines(Data.noisy)
  g=g+1
  
  i=i+step
  
  print(i)
  #lines(Data,col="red",lwd=2)
}

mynewdata_smcabc_modi<-sapply(as.data.frame(saved_estimations_smcabc_modi),function(i) quantile(i, v1))
mynewdata0_smcabc_modi<-mynewdata_smcabc_modi[1,]
mynewdata100_smcabc_modi<-mynewdata_smcabc_modi[7,]

polygon(c(1:length(grnd.truth.I_),length(grnd.truth.I_):1),c(mynewdata0_smcabc_modi,rev(mynewdata100_smcabc_modi)),col=rgb(0, 0, 1,0.5), border=NA)




#plot
plot(grnd.truth.I_, type="l", xlab = "Time(weeks)",
     ylab = "Proportion infected",ylim=c(0,0.36),col="white",cex.axis=1.2)

lines(mynewdata0_mcmc,col="red", lwd=1.5,lty=5)
lines(mynewdata100_mcmc,col="red", lwd=1.5,lty=5)

lines(mynewdata0_smcabc_modi,col="blue",lty=5, lwd=1.5)
lines(mynewdata100_smcabc_modi,col="blue",lty=5, lwd=1.5)

lines(mynewdata0_smcabc_curr,col="green",lty=5, lwd=1.5)
lines(mynewdata100_smcabc_curr,col="green",lty=5, lwd=1.5)



# Finally plotting the posterior mean LOESS line

points(data.rep,col = "black",pch=16)

text=c("True data", "MCMC", "Modified SMC ABC","SMC ABC")
legend("topright", legend = c("True data", "MCMC", "Modified SMC ABC","SMC ABC"),
       col = c("black","red","blue", "green" ),  lty=c(NA,5,5,5),pch=c(16,NA,NA,NA),
       ncol=1,pt.cex=2,border = NA,text.width = strwidth(text)[3])




