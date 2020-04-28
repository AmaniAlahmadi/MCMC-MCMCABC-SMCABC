#scatter map for joint posterior MCMC

library(emdbook)

par(mfrow = c(1,2))

MCMC_data_matrix=matrix(NA,nrow=length(MCMC.run[['gamma']]),ncol=2)
MCMC_data_matrix[,1]=MCMC.run[['gamma']]
MCMC_data_matrix[,2]=MCMC.run[['beta']]

plot(MCMC_data_matrix, col=1,xlab=expression(atop(paste(gamma)," \n MCMC Joint Posterior")),ylab=expression(beta),xlim=c(0.320,0.350),ylim=c(0.87,0.93))
HPDregionplot(MCMC_data_matrix, prob = c(0.95, 0.75, 0.5), col=c("red", "blue", "green"), lwd=3, add=TRUE)
legend("bottomright", legend = c("95%", "75%", "50%"), col = c("red", "blue", "green"), lty=c(1,1,1), lwd=c(3,3,3))


## scatter map for joint posterior for  ABC_SMC
MCMCABC_smc_data_matrix=matrix(NA,nrow=length(ABC_SMC.run[,1]),ncol=2)
MCMCABC_smc_data_matrix[,1]=ABC_SMC.run[,1]
MCMCABC_smc_data_matrix[,2]=ABC_SMC.run[,2]
plot(MCMCABC_smc_data_matrix,xlab=expression(atop(paste(gamma)," \n ABC SMC Joint Posterior")),ylab=expression(beta),xlim=c(0.320,0.350),ylim=c(0.87,0.93))
HPDregionplot(MCMCABC_smc_data_matrix, prob = c(0.95, 0.75, 0.5), col=c("red", "blue", "green"), lwd=3, add=TRUE)
legend("bottomright", legend = c("95%", "75%", "50%"), col = c("red", "blue", "green"), lty=c(1,1,1), lwd=c(3,3,3))

