set.seed(123)
arg=1

#___________________________________________________________

source("ABC_MCMC_SIR_final.R", echo = TRUE)
source("Initial_parameters_final.R", echo = TRUE)

#___________________________________________________________

                          
########################
# run MCMCM on data
########################
start_time <- proc.time()

MCMC.run <- analyze.data_MCMC(data.rep = data.rep,
                                theta.true = theta.true,
                                analyze.data.inputs = analyze.data.inputs)
                               
total_time_McMc <- proc.time() - start_time


# Save the samples and plot the trace plots __________________________________

filename <- paste0("MCMC.run",1,".csv")
write.csv(MCMC.run,filename)

par(mfrow = c(3,1), mar=c(1,1,1,1))
plot(MCMC.run[['gamma']],type="l", main = expression(paste(gamma," trace plot")), 
     ylab = "probability density", xlab = expression(gamma),
     pch = NA,cex = 0.25, pin = c(13, 3.25))
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
abline(h =theta.true['gamma'], col="red" )


plot(MCMC.run[['beta']],type="l", main = expression(paste(beta," trace plot")), 
     ylab = "Accepted Value", xlab = "Iteration number", 
     pch = NA, cex = 0.25, pin = c(13, 3.25))
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
abline(h =theta.true['beta'], col="red" )

plot(MCMC.run[['sigma2']],type="l", main = expression(paste(sigma^2," trace plot")), 
     ylab = "Accepted Value", xlab = "Iteration number", 
     pch = NA, cex = 0.25, pin = c(13, 3.25))
grid(nx = NA,  ny = NULL, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
abline(h =theta.true['sigma2'], col="red" )

title("Trace plots", outer = TRUE)

# ACF plots ____________________________________________________

acf(MCMC.run[['gamma']])
acf(MCMC.run[['beta']])
acf(MCMC.run[['sigma2']])

title("ACF plots", outer = TRUE)


# Convergence curves plots ____________________________________________________

par(mfrow = c(3,1),  mar=c(1,1,1,1))
plot(density(MCMC.run[['beta']]),col=1, main = expression(paste("Marginal posterior of ", beta)),
     ylab = "probability density", xlab = expression(beta),
     xlim = c(beta_true-2*beta.cand.sd, beta_true+2*beta.cand.sd),
     type="l", axes=FALSE)
abline(v =0.9, col="red" )
axis(1, seq(beta_true-2*beta.cand.sd, beta_true+2*beta.cand.sd, beta.cand.sd))

plot(density(MCMC.run[['gamma']]),col=1, main = expression(paste("Marginal posterior of ", gamma)), 
     ylab = "probability density", xlab = expression(gamma),
     xlim = c(gamma_true-2*gamma.cand.sd, gamma_true+2*gamma.cand.sd),
     type="l", axes=FALSE)
abline(v =1/3, col="red" )
axis(1, seq(gamma_true-2*gamma.cand.sd, gamma_true+2*gamma.cand.sd, gamma.cand.sd))##to change

plot(density(MCMC.run[['sigma2']]),col=1, main = expression(paste("Marginal posterior of ", sigma^2)),
     ylab = "Accepted Value", xlab = "Iteration number", 
     xlim = c(sigma2_true-4*sigma2.cand.sd, sigma2_true+4*sigma2.cand.sd),
     type="l", axes=FALSE)
abline(v =(0.0001), col="red" )
axis(1, seq(sigma2_true-4*sigma2.cand.sd, sigma2_true+4*sigma2.cand.sd, sigma2.cand.sd))##to change
title("Convergence curves", outer = TRUE)


#plot the data and fitted-------------------------------
plot(out[,"time"], out[,"I"],type="l", main = "Data and fit", ylab = "Number infection cases", xlab = "Time (weeks)", col = "red",ylim=c(0,0.32),lwd=2)
lines(out[min_time:max_time,"time"],data.rep, type="p",col = "blue",pch = 8)
# fitted ODE solution for infection profile
gamma_median_MCMC = median(MCMC.run[['gamma']])
beta_median_MCMC = median(MCMC.run[['beta']])
parmsTemp_MCMC= c(beta= beta_median_MCMC,gamma= gamma_median_MCMC)
DataTemp_MCMC <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parmsTemp_MCMC)))
lines(DataTemp_MCMC[,"time"], DataTemp_MCMC[,"I"], type="l",col = "green")


#scatter plots pairs
parm.matrix=matrix(NA,nrow=length(MCMC.run[['gamma']]),ncol=3)
parm.matrix[,1]=MCMC.run[['gamma']]
parm.matrix[,2]=MCMC.run[['beta']]
parm.matrix[,3]=MCMC.run[['sigma2']]
colnames(parm.matrix)=c('gamma','beta','sigma2')
pairs(parm.matrix,col="light blue")


########################
# run ABC-SMC on data
########################
start_time <- proc.time()

ABC_SMC.run <- analyze.data_ABC_smc(data.rep = data.rep,
                            theta.true = theta.true,
                            analyze.data.inputs = analyze.data.inputs)

total_time_abc_smc <- proc.time() - start_time


# Save the samples __________________________________

filename <- paste0("ABC_SMC.run_",arg,".csv")
write.csv(ABC_SMC.run,filename)

#scatter plots pairs
parm.matrix=matrix(NA,nrow=length(ABC_SMC.run[,1]),ncol=2)
parm.matrix[,1]=ABC_SMC.run[,1]
parm.matrix[,2]=ABC_SMC.run[,2]
colnames(parm.matrix)=c('gamma','beta')
pairs(parm.matrix,col="light blue")



#the mean absolute errors
# MCMC 
samples_mcmc_beta <- MCMC.run[,2]
samples_mcmc_gamma <- MCMC.run[,1]


samples_ABC_smc_beta_R=ABC_SMC.run[,2]
samples_ABC_smc_gamma_R=ABC_SMC.run[,1]

#Function that returns Mean Absolute Error
MAE_FUNCTION=function(par,length,samples){
  actual_beta <- c(rep(par,length))
  error=actual_beta-samples
  Mean_Absolute_Error=mean(abs(error))
  return(Mean_Absolute_Error)
}

#beta
Mean_Absolute_Error_MCMC_beta=MAE_FUNCTION(0.9,length(samples_mcmc_beta),samples_mcmc_beta)
Mean_Absolute_Error_ABCSMC_R_beta=MAE_FUNCTION(0.9,length(samples_ABC_smc_beta_R),samples_ABC_smc_beta_R)

#gamma
Mean_Absolute_Error_MCMC_gamma=MAE_FUNCTION((1/3),length(samples_mcmc_gamma),samples_mcmc_gamma)
Mean_Absolute_Error_ABCSMC_R_gamma=MAE_FUNCTION((1/3),length(samples_ABC_smc_gamma_R),samples_ABC_smc_gamma_R)


#estimated value:
beta_mcmc_median=median(samples_mcmc_beta)
gamma_mcmc_median=median(samples_mcmc_gamma)

beta_R_median=median(samples_ABC_smc_beta_R)
gamma_R_median=median(samples_ABC_smc_gamma_R)





