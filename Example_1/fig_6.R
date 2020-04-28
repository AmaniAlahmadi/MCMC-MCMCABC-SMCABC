

source("ABC_MCMC_SIR_final.R", echo = TRUE)
source("Initial_parameters_final.R", echo = TRUE)

#generat data with different noise_____________

sigma2_true_1=0.0001
data.rep_1=generat.noisy.data(grnd.truth.I_,sigma2_true_1)


sigma2_true_2=0.0005
data.rep_2=generat.noisy.data(grnd.truth.I_,sigma2_true_2)


sigma2_true_3=0.001
data.rep_3=generat.noisy.data(grnd.truth.I_,sigma2_true_3)



# Run MCMC_______________________________________
MCMC.run_1 <- analyze.data_MCMC(data.rep = data.rep_1,
                              theta.true = theta.true,
                              analyze.data.inputs = analyze.data.inputs)
MCMC.run_2 <- analyze.data_MCMC(data.rep = data.rep_2,
                              theta.true = theta.true,
                              analyze.data.inputs = analyze.data.inputs)
MCMC.run_3 <- analyze.data_MCMC(data.rep = data.rep_3,
                              theta.true = theta.true,
                              analyze.data.inputs = analyze.data.inputs)

#_________________________________________________
dens_gamma_1=density(MCMC.run_1[,'gamma'])
dens_gamma_2=density(MCMC.run_2[,'gamma'])
dens_gamma_3=density(MCMC.run_3[,'gamma'])
dens_beta_1=density(MCMC.run_1[,'beta'])
dens_beta_2=density(MCMC.run_2[,'beta'])
dens_beta_3=density(MCMC.run_3[,'beta'])

plot(dens_gamma_3)
lines(dens_gamma_1)
lines(dens_gamma_2)

#Run ABC SMC______________________________________

ABC_SMC.run_1 <- analyze.data_ABC_smc(data.rep = data.rep_1,
                                    theta.true = theta.true,
                                    analyze.data.inputs = analyze.data.inputs)

ABC_SMC.run_2 <- analyze.data_ABC_smc(data.rep = data.rep_2,
                                    theta.true = theta.true,
                                    analyze.data.inputs = analyze.data.inputs)
ABC_SMC.run_3 <- analyze.data_ABC_smc(data.rep = data.rep_3,
                                    theta.true = theta.true,
                                    analyze.data.inputs = analyze.data.inputs)


#_______________________________________________

densabc_gamma_1=density(ABC_SMC.run_1[,1])
densabc_gamma_2=density(ABC_SMC.run_2[,1])
densabc_gamma_3=density(ABC_SMC.run_3[,1])
densabc_beta_1=density(ABC_SMC.run_1[,2])
densabc_beta_2=density(ABC_SMC.run_2[,2])
densabc_beta_3=density(ABC_SMC.run_3[,2])


#Plot the results________________________________________________

par(mfrow = c(2,3))
#gamma___________________________________________________________
plot(dens_gamma_1,ylim=c(0,2000),xlim=c(0.29,0.37),col="white",
     xlab=expression(gamma),main=expression(paste(sigma^2,"=0.0001")))

polygon(densabc_gamma_1, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_gamma_1, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.gamma, col=1 )

plot(dens_gamma_2,ylim=c(0,2000),xlim=c(0.29,0.37),
     col="white",xlab=expression(gamma),main=expression(paste(sigma^2,"=0.0005")))
polygon(densabc_gamma_2, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_gamma_2, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.gamma, col=1 )

plot(dens_gamma_3,ylim=c(0,2000),xlim=c(0.29,0.37),
     col="white",xlab=expression(gamma),main=expression(paste(sigma^2,"=0.001")))
polygon(densabc_gamma_3, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_gamma_3, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.gamma, col=1 )
#beta________________________________________________________________

plot(densabc_beta_1,ylim=c(5,2000),xlim=c(0.87,0.93),
     col="white",xlab=expression(beta),main=NA)
polygon(densabc_beta_1, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_beta_1, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.beta, col=1 )

plot(densabc_beta_2,ylim=c(5,2000),xlim=c(0.87,0.93),
     col="white",xlab=expression(beta),main=NA)
polygon(densabc_beta_2, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_beta_2, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.beta, col=1 )

plot(dens_beta_3,ylim=c(5,2000),xlim=c(0.87,0.93),
     col="white",xlab=expression(beta),main=NA)
polygon(densabc_beta_3, col=rgb(0.4, 0.9, 0,0.7), border=NA)
polygon(dens_beta_3, col=rgb(1, 0, 0,0.7), border=NA)
abline(v =initial.beta, col=1 )












