
source("ABC_MCMC_SIR_final.R", echo = TRUE)
source("Initial_parameters_final.R", echo = TRUE)

#Plot fig 10 on the paper_______________________


Vaart_method=analyze.data_Vaart_method(data.rep,theta.true,
                                       analyze.data.inputs)


MCMC.run<- analyze.data_MCMC(data.rep = data.rep,
                             theta.true = theta.true,
                             analyze.data.inputs = analyze.data.inputs)


plot(Vaart_method[,'gamma.samples'],Vaart_method[,'beta.samples'],
     col=2,xlab=expression(gamma),ylab=expression(beta))
points (MCMC.run[,'gamma'],MCMC.run[,'beta'])