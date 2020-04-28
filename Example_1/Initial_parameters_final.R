# Load libraries
library(deSolve)
library(rootSolve)
library(mvtnorm )
library(MCMCpack)
library(MASS)

#initial.param.inputs______________________________________
initial.beta <-0.9
initial.gamma <- 1/3
initial.sigma2 <- 0.0001


initial.param.inputs <- c(initial.beta = initial.beta,
                          initial.gamma = initial.gamma,
                          initial.sigma2 = initial.sigma2)

#generate.data.inputs_______________________________________
T1=50                                       #real time[0,50]
del=0.5                                          #integration interval
true_time=seq(0,T1,by= del)                 #true times

N=1
I=1.27*10^-6
R=0
S=N-I-R
ICS= c(S=S,I=I,R=R)                         #init vals of latent states

#generate.data______________________________________________
##Ground Truth
ODE_system<-function(true_time, ICS, theta.true) {
  with(as.list(c(ICS, theta.true)),{
    dS=-beta*S*I
    dI=beta*S*I-gamma*I
    dR=gamma*I
    
    der <-c(dS,dI,dR)
    # return the rate of change
    return(list(der))
  })
}
beta_true=initial.param.inputs[['initial.beta']]
gamma_true=initial.param.inputs[['initial.gamma']]
sigma2_true=initial.param.inputs[['initial.sigma2']]
theta.true= c(beta=beta_true,gamma=gamma_true,sigma2= sigma2_true)
out <-ode(y = ICS, times = true_time, func = ODE_system, parms = theta.true)


#save ground truth I for window
grnd.truth.I <- out[,3]
min_time <- min(which((grnd.truth.I>0.01)==TRUE))
max_time <- max(which((grnd.truth.I>0.01)==TRUE))
if (min_time==Inf){
  min_time <- min(which((grnd.truth.I>0.001)==TRUE))
  max_time <- max(which((grnd.truth.I>0.001)==TRUE))
}
if (min_time==Inf){
  min_time <- 0
  max_time <- length(grnd.truth.I)
}

grnd.truth.I_ <- grnd.truth.I[min_time:max_time]

generate.data.inputs <- c(sigma2_true,ICS,min_time,max_time,true_time)

#add noise, but make sure it is not negative
generat.noisy.data=function(grnd.truth.I_,sigma2_true){
for(i in 1:100000){
  data.rep <- grnd.truth.I_ +
    rnorm(length(grnd.truth.I_),0,sqrt(sigma2_true))
  if (sum(data.rep<0)==0){
    break
  }
}
  return(data.rep)
}
data.rep=generat.noisy.data(grnd.truth.I_,sigma2_true)
#Save the data
#saveRDS(data.rep,"Sir_data.rds")
#read the data
data.rep=readRDS("Sir_data.rds")
#discrepancy function_______________________________________
rho<-function(x, y){sqrt(sum((x-y)^2))}

rho.new=function(x, y,sigma.value){sqrt(sum(((x-y)/sigma.value)^2))}

#analyze.data.inputs________________________________________

#Candidate search standard deviations
beta.cand.sd <- 0.04
gamma.cand.sd <- 0.04
sigma2.cand.sd <- 0.0001

Iterations <- 13000
burn_in <- 5000
rounds=11
number.of.particals=20000  #for Vaart method change t0 2*10^6

#Define Prior functions
ODE_log.Prior_gamma=function(x){dunif(x,0,2,log=T)}
ODE_log.Prior_beta=function(x){dunif(x,0,2,log=T)}
ODE_log.Prior_sigma2=function(x){log(1/dgamma(x,1,1))}


#Concatenate analyze.data.inputs
analyze.data.inputs <- list(beta.cand.sd=beta.cand.sd,
                            gamma.cand.sd=gamma.cand.sd,
                            sigma2.cand.sd=sigma2.cand.sd,
                            ICS=ICS,
                            true_time=true_time,
                            min_time=min_time,
                            max_time=max_time,
                            iterations=Iterations,
                            burn_in=burn_in,
                            rounds=rounds,
                            number.of.particals=number.of.particals)


