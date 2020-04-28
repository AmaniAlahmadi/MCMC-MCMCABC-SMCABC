


source("ABC_MCMC_SIR_final.R", echo = TRUE)
source("Initial_parameters_final.R", echo = TRUE)

#### Figure 5 ####


ABC_PMC=analyze.data_ABC_smc_adaptive_distance(data.rep,theta.true,
                                               analyze.data.inputs)


dat <- matrix(NA,nrow=length(ABC_PMC[1,,1]), ncol = 2)
dat[,1]=ABC_PMC[2,,15]
dat[,2]=ABC_PMC[1,,15]
ch <- chull(dat)
coords <- dat[c(ch, ch[1]), ]  # closed polygon

par(mgp=c(2.3,1,0))# the position of the plot

plot(dat,col=NA,xlab=expression(gamma),ylab=expression(beta),cex.lab=1.7,xlim=c(0.29,0.37),ylim=c(0.85,0.94))

lines(coords, col="red",lwd=2)

points(MCMC.run[['gamma']],MCMC.run[['beta']],col=6)


dat <- matrix(NA,nrow=length(ABC_PMC[2,,6]), ncol = 2)
dat[,1]=ABC_PMC[2,,19]
dat[,2]=ABC_PMC[1,,19]
ch <- chull(dat)
coords <- dat[c(ch, ch[1]), ]  # closed polygon

lines(coords, col=1,lwd=2)



dat <- matrix(NA,nrow=length(ABC_PMC[2,,30]), ncol = 2)
dat[,1]=ABC_PMC[2,,30]
dat[,2]=ABC_PMC[1,,30]
ch <- chull(dat)
coords <- dat[c(ch, ch[1]), ]  # closed polygon



lines(coords, col="blue",lwd=2)
abline(v=0.3333,lty=2,col=5)
abline(h=0.9,lty=2,col=5)

legend("bottomright", legend = c("Round 6", "Round 16","Round 30","MCMC posterior",
                                 "Parameters true values"), col = c("red", "black","blue",6,5),
       lty=c(1,1,1,NA,2), pch=c(NA,NA,NA,1,NA), lwd=c(2,2,2,2,2),cex = 0.9)




