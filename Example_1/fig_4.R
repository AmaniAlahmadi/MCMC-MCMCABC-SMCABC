# Plotting a scatter of posterior lines to represent the uncertainty.
#ploting the uncertinity of the posterior
# Plotting the data
par(mfrow = c(1,1))

plot(data.rep, type="l", xlab = "Time(weeks)",
     ylab = "Proportion infected",ylim=c(0,0.35),xaxt="n")
axis(side=1,at=c(0,20,40,60,80,100),labels=c("0","10","20","30","40","50"))

step=7 
saved_estimations_mcmc=matrix(data=NA,ncol=length(data.rep),nrow=length(MCMC.run[['gamma']])/step)

#MCMC
i=1
g=0
while(i < length(MCMC.run[['gamma']])-4) {
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
mynewdata_mcmc<-sapply(as.data.frame(saved_estimations_mcmc),FUN=quantile)
mynewdata0_mcmc<-mynewdata_mcmc[1,]
mynewdata100_mcmc<-mynewdata_mcmc[5,]

polygon(c(1:length(data.rep),length(data.rep):1),c(mynewdata0_mcmc,rev(mynewdata100_mcmc)),col="pink",border=NA)



#__________________________________________________________________
#ABC_SMC
i=1
g=0
step=1 
saved_estimations_abc_smc=matrix(data=NA,ncol=length(data.rep),nrow=length(ABC_SMC.run[,1])/step)

while(i <= length(ABC_SMC.run[,1])) {
  gamma_pre=(ABC_SMC.run[,1][i])
  beta_pre=(ABC_SMC.run[,2][i])
  parms= c(beta= beta_pre,gamma=gamma_pre)
  Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
  Data= Data[,3]
  Data= Data[min_time:max_time]
  
  saved_estimations_abc_smc[g+1,]=Data
  g=g+1
  i=i+step
  
  print(i)
  lines(Data,col="red",lwd=2)
}


mynewdata_abc_smc<-sapply(as.data.frame(saved_estimations_abc_smc),FUN=quantile)
mynewdata0_abc_smc<-mynewdata_abc_smc[1,]
mynewdata100_abc_smc<-mynewdata_abc_smc[5,]

polygon(c(1:length(data.rep),length(data.rep):1),c(mynewdata0_abc_smc,rev(mynewdata100_abc_smc)),col="red",border=NA)



# Finally plotting the posterior mean LOESS line
lines(data.rep,type="l", col = "blue", xlab = "....",
      ylab = "....", main = "....")
x=seq(0,100,by=1)
points(data.rep,col = "blue",pch=16)


legend("topleft", legend = c("True data", "MCMC", "SMC ABC"),
       col = c("blue", "pink", "red"),  pch=c(16,15,15),ncol=3,pt.cex=3)




