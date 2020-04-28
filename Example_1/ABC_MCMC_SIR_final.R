

#analyze.data_______________________________________________
analyze.data_MCMC <- function(data.rep,theta.true,
                              analyze.data.inputs){
  
  #Decompose analyze.data.inputs
  beta.cand.sd <- analyze.data.inputs[['beta.cand.sd']]
  gamma.cand.sd <- analyze.data.inputs[['gamma.cand.sd']]
  sigma2.cand.sd <- analyze.data.inputs[['sigma2.cand.sd']]
  ICS <- analyze.data.inputs[['ICS']]
  true_time <- analyze.data.inputs[['true_time']]
  min_time <- analyze.data.inputs[['min_time']]
  max_time <- analyze.data.inputs[['max_time']]
  iterations <- analyze.data.inputs[['iterations']]
  burn_in <- analyze.data.inputs[['burn_in']]
  
 #get good starting guess for beta and gamma
  temp_index = which(data.rep ==max(data.rep))
  Imax = data.rep[temp_index]
  tmax = true_time[seq(min_time,max_time)[temp_index]]
  
  mod = lm(log(data.rep[temp_index:length(data.rep)]) ~ true_time[seq(min_time,max_time)[temp_index]:max_time])
  gamma_0 = -coef(mod)[[2]]
  Imax_fun <- function (val) Imax + val - val*log(val) -1
  inv_ratio <- uniroot(Imax_fun, c(0.000001, 1))$root
  ratio = 1/inv_ratio
  beta_0 =ratio*gamma_0
  sigma2_0=theta.true[['sigma2']]  
  
  ## old way to get initial guesses
  #beta_0=theta.true[['beta']]
  #gamma_0=theta.true[['gamma']]
  #sigma2_0=theta.true[['sigma2']]
  
  output.df<- data.frame(gamma=NA, beta=NA, sigma2=NA, cand.ratio=NA, cand.like=NA, cand.prior=NA)
  
  pb <- txtProgressBar(min = 0, max = iterations, initial = 1,
  style =3 ,
  char = "=")
  
  parms= c(beta= beta_0,gamma= gamma_0)
  Data <-ode(y = ICS, times = true_time, func = ODE_system, parms = parms)
  Data <- Data[,3]
  Data <- Data[min_time:max_time]
  
  curnt.prior.log = ODE_log.Prior_gamma(gamma_0)+ODE_log.Prior_beta(beta_0)+ODE_log.Prior_sigma2(sigma2_0)
  curnt.like.log = sum(dnorm(data.rep,mean=Data, sqrt(sigma2_0),log=T))
  
  
  #counting the acceptance value	
  countf_gamma_beta=0
  countf_sigma2=0
  x=0      #count of itreations
  d=0      #index
  #Start the loop
  while(countf_gamma_beta<iterations){
    
    setTxtProgressBar(pb,countf_gamma_beta+1)
    
    ###########################################
    #   Begin estimation for gamma and beta   #
    ###########################################
    
    repeat {
      gamma.cand<-gamma_0+rnorm( n=1,mean=0,sd=gamma.cand.sd) 
      beta.cand<-beta_0 +rnorm( n=1,mean=0,sd=beta.cand.sd) 
      
        if ( (beta.cand>gamma.cand) ){
        break
    }
    }
    cand.prior.log = ODE_log.Prior_gamma(gamma.cand)+ODE_log.Prior_beta(beta.cand)+ODE_log.Prior_sigma2(sigma2_0)
    parms= c(beta= beta.cand,gamma=gamma.cand)
    Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
    Data= Data[,3]
    Data= Data[min_time:max_time]
    
    cand.like.log =sum(dnorm(data.rep,mean=Data, sqrt(sigma2_0),log=T))
    
    #calculate the ratio   
    cand.ratio<- (cand.prior.log+cand.like.log-
                    curnt.prior.log-curnt.like.log)
    
    if (log(runif(1))< cand.ratio){
      #store the accepted values
      gamma_0<-gamma.cand
      beta_0<-beta.cand
      
      curnt.prior.log<-cand.prior.log
      curnt.like.log<-cand.like.log
      
      countf_gamma_beta = countf_gamma_beta +1
    }
    
    ###################################
    #   Begin estimation for sigma2   #
    ###################################
    
    repeat {
      sigma2.cand<-sigma2_0 +rnorm(n=1,mean=0,sd=sigma2.cand.sd) 
      
      if (sigma2.cand>0)  #To make sure sigma2.cand>0
        break
    }
    
    cand.prior.log = ODE_log.Prior_gamma(gamma_0)+ODE_log.Prior_beta(beta_0)+ODE_log.Prior_sigma2(sigma2_0)
    
    parms= c(beta= beta_0,gamma=gamma_0)
    Data <-((ode(y = ICS, times = true_time, func = ODE_system, parms = parms)))
    Data= Data[,3]
    Data= Data[min_time:max_time]
    
    cand.like.log =sum(dnorm(data.rep,mean=Data, sqrt(sigma2.cand),log=T))
    
    #calculate the ratio   
    cand.ratio<- (cand.prior.log+cand.like.log-
                    curnt.prior.log-curnt.like.log)
    if (log(runif(1))< cand.ratio){
      #store the accepted values
      sigma2_0 <- sigma2.cand
      
      curnt.prior.log<-cand.prior.log
      curnt.like.log<-cand.like.log
      
      countf_sigma2 = countf_sigma2 +1
      
    }
    x=x+1
    if (x>500&&x<3000){
      d=d+1
      if((countf_gamma_beta/x)>0.5&d==100){
        d=0
        gamma.cand.sd=gamma.cand.sd*1.2
        beta.cand.sd= beta.cand.sd*1.2
      }
      else if((countf_gamma_beta/x)<0.23&d==100){
        d=0

        gamma.cand.sd=gamma.cand.sd/1.2
        beta.cand.sd= beta.cand.sd/1.2
      }
      
    }
    if (x>burn_in) {
      output.df[x-burn_in,]<-c(gamma_0, beta_0, sigma2_0, cand.ratio, exp(cand.like.log), exp(cand.prior.log))
    }
  }
  print("Number of values accepted")
  print(c("Gamma, Beta: ", countf_gamma_beta))
  print(c("Sigma^2: ", countf_sigma2))
  print(c("Number of itreations: ", x))
  
  return(output.df)
  
}




#___________________________________________________________

analyze.data_ABC_smc <- function(data.rep,theta.true,
                                 analyze.data.inputs){
  
  ICS <- analyze.data.inputs[['ICS']]
  true_time <- analyze.data.inputs[['true_time']]
  min_time <- analyze.data.inputs[['min_time']]
  max_time <- analyze.data.inputs[['max_time']]
  rounds<- analyze.data.inputs[['rounds']]
  number.of.particals<- analyze.data.inputs[['number.of.particals']]
  
  
  
  
  ####
  # Number of SMC rounds
  rounds <- rounds
  
  # Number of iterations per round
  N <- number.of.particals
  
  # Fraction of samples to save per round
  tau <- c(0.1,rep(0.25,rounds-1))
  
  
  # Number of samples to save from first round
  Q <- N*tau[1]
  
  # Threshold value for simulated vectors to be saved
  epsilon <-rep(NA,rounds)
  # Storage vector for dist measure from simulated data
  dist_vals <- rep(NA,N)
  # Draw initial beta and gamma  samples for first round from prior
  j=1
  beta_sam=c(NA)
  gamma_sam=c(NA)

  while(j <= N){
    generate.beta_ <-runif(1,0,2)
    generate.gamma_ <-runif(1,0,2)

    if((generate.beta_>generate.gamma_)){
      beta_sam[j]=generate.beta_
      gamma_sam[j]=generate.gamma_
      j=j+1
    }
  }

  
  
  #counting the acceptance value	
  countf_gamma_beta=0
  
  ### -------- 1st round -------- ##
  print("Starting the first stage")
  pb <- txtProgressBar(min = 0, max = N, initial = 1,
                       style =3 ,
                       char = "=")
  # Draw N data sets and compare
  for(i in 1:N){
    parms= c(beta= beta_sam[i],gamma=gamma_sam[i])
    #Solve the ODE       
    Data <-(ode(y = ICS, times = true_time, func = ODE_system, parms = parms))
    Data= Data[,3]
    Data= Data[min_time:max_time]
    dist_vals[i] <- rho(Data,data.rep)
    
    setTxtProgressBar(pb,i+1)
    
  }
  # Sort dist values and save the top Q beta and gamma
  dist_indexes <- sort(dist_vals, index.return=T)
  save_indexes <- dist_indexes$ix[1:Q]
  saved_beta <- beta_sam[save_indexes]
  saved_gamma <- gamma_sam[save_indexes]
  plot(saved_gamma,saved_beta,col=1,xlab="gamma",ylab="beta")
  #initial value for kernal width
  beta.sd<-1/2*(max(saved_beta)-min(saved_beta))
  gamma.sd<-1/2*(max(saved_gamma)-min(saved_gamma))
  # first epsilon is the max dist value
  epsilon[1] <- dist_indexes$x[Q] 
  
  #initial wieghts for all paramters
  w_beta=rep(1,length(saved_beta))/length(saved_beta)
  w_gamma=rep(1,length(saved_beta))/length(saved_beta)
 
  ## -------- 2nd through rth rounds -------- ##
  print("Starting the second stage")
  pb <- txtProgressBar(min = 0, max = rounds, initial = 1,
                       style =3 ,
                       char = "=")
  for(r in 2:rounds){
    curr_num_saved <- 0
    dist_vals <- rep(NA,Q)
    dist_vals_sigma<- rep(NA,Q)
    tmp_saved_beta <- rep(NA,Q)
    tmp_saved_gamma <- rep(NA,Q)
    tmp_weights_beta<-rep(NA,Q)
    tmp_weights_gamma<-rep(NA,Q)

    while(curr_num_saved < Q){
      
      repeat{
        curr_beta <- sample(saved_beta,1,prob=w_beta)
        curr_beta <- curr_beta + rnorm(1,0,beta.sd)
        curr_gamma <- sample(saved_gamma,1,prob=w_gamma)
        curr_gamma <- curr_gamma + rnorm(1,0,gamma.sd)
        if(curr_beta>curr_gamma){break}
      }
      
      
      
      cand.prior.log<-ODE_log.Prior_gamma(curr_gamma)+ODE_log.Prior_beta(curr_beta)
      if (exp(cand.prior.log)>0 ){
        parms= c(beta= curr_beta,gamma=curr_gamma)
        
        #Solve the ODE  
        Data <-(ode(y = ICS, times = true_time, func = ODE_system, parms = parms))
        Data= Data[,3]
        curr_dat= Data[min_time:max_time]
       
        if (sum(curr_dat<0)==0) {

          curr_dist <- rho(data.rep,curr_dat)
          
          if(curr_dist < epsilon[r-1]){
            curr_num_saved <- curr_num_saved + 1
            dist_vals[curr_num_saved] <- curr_dist
            tmp_saved_beta[curr_num_saved] <- curr_beta
            tmp_saved_gamma[curr_num_saved] <- curr_gamma
            #find the wight for beta
            weights.numerator_beta<- exp(ODE_log.Prior_beta(curr_beta))
            weights.denominator_beta<- sum(w_beta*dnorm(curr_beta,saved_beta,beta.sd))
            tmp_weights_beta[curr_num_saved]<-weights.numerator_beta/weights.denominator_beta
            #find the wight for gamma
            weights.numerator_gamma<- exp(ODE_log.Prior_gamma(curr_gamma))
            weights.denominator_gamma<- sum(w_gamma*dnorm(curr_gamma,saved_gamma,gamma.sd))
            tmp_weights_gamma[curr_num_saved]<-weights.numerator_gamma/weights.denominator_gamma
          
            }
        }
      }
    }
    setTxtProgressBar(pb,r+1)
    
    # Save beta and gamma samples
    dist_indexes <- sort(dist_vals, index.return=T)
    save_indexes <- dist_indexes$ix[1:(tau[r]*Q)]
    

    epsilon[r] <- dist_indexes $x[(tau[r]*Q)]
    saved_beta <- tmp_saved_beta[save_indexes]
    saved_gamma<- tmp_saved_gamma[save_indexes]
    distanc=dist_vals[save_indexes]
    
    
    
    beta.sd.number = 1/2*(max(saved_beta)-min(saved_beta))
    gamma.sd.number = 1/2*(max(saved_gamma)-min(saved_gamma))

    beta.sd= beta.sd.number
    gamma.sd= gamma.sd.number
    
    w_beta=tmp_weights_beta[save_indexes]
    w_beta=w_beta/sum(w_beta)
    w_gamma=tmp_weights_gamma[save_indexes]
    w_gamma=w_gamma/sum(w_gamma)
    points(saved_gamma,saved_beta,col=r)
    
    
  }
  
  #store final samples of parameters
  output.df<- matrix(NA,nrow=length(saved_beta), ncol=2)
  
  output.df[,1]<-c(saved_gamma)
  output.df[,2]<-c(saved_beta)

  
  print("Number of values accepted")
  print(c("Gamma, Beta: ", countf_gamma_beta))
  print(epsilon)
  return(output.df)
  
}




#ABC_SMC with adaptive_distance_____________________________
analyze.data_ABC_smc_adaptive_distance <- function(data.rep,theta.true,
                                 analyze.data.inputs){
  ICS <- analyze.data.inputs[['ICS']]
  true_time <- analyze.data.inputs[['true_time']]
  min_time <- analyze.data.inputs[['min_time']]
  max_time <- analyze.data.inputs[['max_time']]
  rounds<-30# analyze.data.inputs[['rounds']]
  number.of.particals<- analyze.data.inputs[['number.of.particals']]
  
  
  
  
  ####
  # Number of SMC rounds
  rounds <- rounds
  
  # Number of iterations per round
  N <- number.of.particals
  
  # Fraction of samples to save per round
  tau <- 0.2
  
  
  # Number of samples to save from first round
  Q <- N*tau[1]
  
  
  # Threshold value for simulated vectors to be saved
  epsilon <- rep(NA,rounds)
  # Storage vector for dist measure from simulated data
  dist_vals <- rep(NA,N)
  # Draw initial gamma and beta samples for first round from prior
  j=1
  beta_sam=c(NA)
  gamma_sam=c(NA)
  error.valus=c(NA)
  while(j <= N){
    generate.beta_v <-runif(1,0,2)
    generate.gamma_v <-runif(1,0,2)
    beta_sam[j]=generate.beta_v
    gamma_sam[j]=generate.gamma_v
    j=j+1 
  }
  
  
  #counting the acceptance value	
  countf_gamma_beta=0
  
  ### -------- 1st round -------- ##
  print("Starting the first stage")
  pb <- txtProgressBar(min = 0, max = N, initial = 1,
                       style =3 ,
                       char = "=")
  #matrix to save all simulation
  simulation.matrix=matrix(NA,nrow=length(data.rep),ncol=N)
  #a  array to save all the itreations
  array.beta.gamma=array(NA, dim = c(2, (tau*Q), rounds))
  sigma.store=c(NA)
  # Draw N data sets and calculate the distance
  for(i in 1:N){
    parms= c(beta= beta_sam[i],gamma=gamma_sam[i])
    #Solve the ODE       
    Data <-(ode(y = ICS, times = true_time, func = ODE_system, parms = parms))
    Data= Data[,3]
    Data= Data[min_time:max_time]
    dist_vals[i] <- rho(Data,data.rep)
    simulation.matrix[,i]=Data
    sigma.store[i]=mad(Data)
    setTxtProgressBar(pb,i+1)
  }
  
  # calculate the new distance_____________________
  new.dist_vals=rep(NA,N)
  for(ii in 1:(N)){
    sigma.value=(sigma.store[ii])
    new.dist_vals[ii]=rho.new(data.rep,simulation.matrix[,ii],sigma.value)
  }
  new.dist_indexes <- sort(new.dist_vals, index.return=T)
  save_indexes <- new.dist_indexes$ix[1:Q]
  quntile=(tau*Q) 
  epsilon[1] <- new.dist_indexes$x[quntile] 
  saved_beta <- beta_sam[save_indexes]
  saved_gamma <- gamma_sam[save_indexes]
  array.beta.gamma[1,,1]=saved_beta[1:(tau*Q)]
  array.beta.gamma[2,,1]=saved_gamma[1:(tau*Q)]
  beta.sd=1/2*(max(saved_beta)-min(saved_beta))
  gamma.sd=1/2*(max(saved_gamma)-min(saved_gamma))
  w_beta=rep(1,length(save_indexes))/length(save_indexes)
  w_gamma=rep(1,length(save_indexes))/length(save_indexes)
  dist_vals_old=dist_vals
  plot(saved_gamma,saved_beta)
  
  ## -------- 2nd through rth rounds -------- ##
  print("Starting the second stage")
  pb <- txtProgressBar(min = 0, max = rounds, initial = 1,
                       style =3 ,
                       char = "=")
  for(r in 2:rounds){
    s=1
    curr_num_saved <- 0
    dist_vals <- rep(NA,Q)
    tmp_saved_beta <- c()
    tmp_saved_gamma <-c()
    tmp_weights_beta<-c()
    tmp_weights_gamma<-c()
    sigma.store=c()
    while(curr_num_saved < Q){
      curr_beta <- sample(saved_beta,1,prob=w_beta)
      curr_beta <- curr_beta + runif(1,-beta.sd,beta.sd)
      curr_gamma <- sample(saved_gamma,1,prob=)
      curr_gamma <- curr_gamma + runif(1,-gamma.sd,gamma.sd)
      
      
      cand.prior.log<-ODE_log.Prior_gamma(curr_gamma)+ODE_log.Prior_beta(curr_beta)
      if (exp(cand.prior.log)>0 ){
        parms= c(beta= curr_beta,gamma=curr_gamma)
        
        #Solve the ODE  
        Data <-(ode(y = ICS, times = true_time, func = ODE_system, parms = parms))
        Data= Data[,3]
        curr_dat= Data[min_time:max_time]
        
        #find the distance 1
        curr_dist <- rho(data.rep,curr_dat)
        #store all the simulations
        simulation.matrix[,s]=curr_dat
        #calculat the scale(mad) for each simulation
        sigma.store[s]=mad(curr_dat)
        
        #Store all the paramters
        tmp_saved_beta[s] <- curr_beta
        tmp_saved_gamma[s] <- curr_gamma
        #store the wight for beta
        weights.numerator_beta<- exp(ODE_log.Prior_beta(curr_beta))
        weights.denominator_beta<- sum(w_beta*dnorm(curr_beta,saved_beta,beta.sd))
        tmp_weights_beta[s]<-weights.numerator_beta/weights.denominator_beta
        #store the wight for gamma
        weights.numerator_gamma<- exp(ODE_log.Prior_gamma(curr_gamma))
        weights.denominator_gamma<- sum(w_gamma*dnorm(curr_gamma,saved_gamma,gamma.sd))
        tmp_weights_gamma[s]<-weights.numerator_gamma/weights.denominator_gamma
        s=s+1
        #accepte the simulation that pass the acceptance rule of the previous iteration
        if(curr_dist < epsilon[r-1]){
          
          curr_num_saved <- curr_num_saved + 1
          dist_vals[curr_num_saved] <- curr_dist
        } 
      }
    }
    setTxtProgressBar(pb,r+1)
    
    #update the distanc (distance2)  
    new.dist_vals=rep(NA,s-1)
    for(ii in 1:(s-1)){
      sigma.value=(sigma.store[ii])
      new.dist_vals[ii]=rho.new(data.rep,simulation.matrix[,ii],sigma.value)
    }
    
    new.dist_indexes <- sort(new.dist_vals, index.return=T)
    save_indexes <- new.dist_indexes$ix[1:(tau*Q)]
    quntile=(tau*Q)
    epsilon[r] <- new.dist_indexes$x[quntile]
    saved_beta <- tmp_saved_beta[save_indexes]
    saved_gamma<- tmp_saved_gamma[save_indexes]
    array.beta.gamma[1,,r]=saved_beta
    array.beta.gamma[2,,r]=saved_gamma
    beta.sd=1/2*(max(saved_beta)-min(saved_beta))
    gamma.sd=1/2*(max(saved_gamma)-min(saved_gamma))
    w_beta=tmp_weights_beta[save_indexes]
    w_beta=w_beta/sum(w_beta)
    w_gamma=tmp_weights_gamma[save_indexes]
    w_gamma=w_gamma/sum(w_gamma)
    dist_vals_old=dist_vals
    points(saved_gamma,saved_beta,col=r)
  }
  
  output.df<- array.beta.gamma
  print("Number of values accepted")
  print(c("Gamma, Beta: ", countf_gamma_beta))
  print(c("epsilon=",epsilon))
  return(output.df)
  
}





#Vaart_method
analyze.data_Vaart_method <- function(data.rep,theta.true,
                                                   analyze.data.inputs){
  ICS <- analyze.data.inputs[['ICS']]
  true_time <- analyze.data.inputs[['true_time']]
  number.of.particals<- analyze.data.inputs[['number.of.particals']]
  
  
  
  
  ###
  # Number of particals
  N <- number.of.particals
  
  # Storage vector for dist measure from simulated data
  dist_vals <- rep(NA,N)
  error.estimate_temp<- rep(NA,N)
  # Draw initial k and din samples for first round from prior
  j=1
  beta_sam=c(NA)
  gamma_sam=c(NA)
  error.valus=c(NA)
  while(j <= N){
    generate.beta_v <-runif(1,0,2)
    generate.gamma_v <-runif(1,0,2)   
   
      beta_sam[j]=generate.beta_v
      gamma_sam[j]=generate.gamma_v
      j=j+1
  }

  
  ### -------- 1st round -------- ##
  simulation.matrix=matrix(NA,nrow=length(data.rep),ncol=N)
  noisy.simulation.matrix=matrix(NA,nrow=length(data.rep),ncol=N)
  
  sigma.store=c(NA)
  store.error=matrix(NA,nrow=length(data.rep),ncol=N)
  print("Starting sampling")
  pb <- txtProgressBar(min = 0, max = N, initial = 1,
                       style =3 ,
                       char = "=")
  # Draw N data sets and calculate the distance to find the best fit and then find the sd
  for(i in 1:N){
    parms= c(beta= beta_sam[i],gamma=gamma_sam[i])
    #Solve the ODE       
    Data <-(ode(y = ICS, times = true_time, func = ODE_system, parms = parms))
    Data= Data[,3]
    Data= Data[min_time:max_time]
    dist_vals[i] <- rho(Data,data.rep)
    simulation.matrix[,i]=Data
    setTxtProgressBar(pb,i+1)
    
  }
  
  # Sort dist values 
  dist_indexes <- sort(dist_vals, index.return=T)

  #find the best fit that minmize the distince bwt data and simulation rho(data.rep,Data)
  best.fit=dist_indexes$ix[1]
  
  #plot the data with the best fit 
  plot(data.rep,col=2,type="l")
  lines(simulation.matrix[,best.fit],col=6)
  
  #find the differnces between the best fit and the data.rep
  e.term=(data.rep-simulation.matrix[,best.fit])
  #find the variance of the first error distrbution
  var.e.term=var(e.term)
  sd.e.term=sqrt(var.e.term)
  
  noisy.observation= data.rep
  
  # calculat the error dinsity
  
  S=c(NA)
  
  for(iii in 1:N){
    
    error.term=simulation.matrix[,iii]-noisy.observation
    
    
    S[iii]=sum((error.term /sd.e.term)^2)
    
  }
  
  
  accept.probability=c(NA)
  
  for(iii in 1:N){
    accept.probability[iii]= dchisq(S[iii],length(x) )*sqrt(S[iii])
  }
  
  c=max(accept.probability) 
  
  normlize.accept.probability=c(NA)
  
  for(iii in 1:N){
    normlize.accept.probability[iii]= accept.probability[iii]/c
  }  
  
  Accepted.normlize.accept.probability=c(NA)
  saved_beta_vec=c(NA)
  saved_gamma_vec=c(NA)
  ss=1
  Accepted.simulation.matrix=matrix(NA,nrow=length(data.rep),ncol=N)
  
  for(iii in 1:N){
    
    if(normlize.accept.probability[iii]>=runif(1)){
      
      Accepted.normlize.accept.probability[ss]= normlize.accept.probability[iii]
      
      saved_beta_vec[ss] <- beta_sam[iii]
      saved_gamma_vec [ss]<- gamma_sam[iii]
      
      Accepted.simulation.matrix[,ss]=simulation.matrix[,iii]
      ss=ss+1
      
      #print(ss)
      
    }
  }
  
  points(saved_gamma_vec,saved_beta_vec)
  #New coverage algorithm for errorcalibrated ABC
  # add nois to all simulation
  dist_vals_noise=c(NA)
  noisy.simulation.matrix=matrix(NA,nrow=length(noisy.observation),ncol=length(saved_beta_vec))
  
  for(k in 1:(ss-1)){
    
    
    noisy.simulation.matrix[,k]=Accepted.simulation.matrix[,k]+rnorm(length(noisy.observation),0,sd.e.term)
    
    dist_vals_noise[k]=rho(noisy.observation, noisy.simulation.matrix[,k])
  }
  
  #find the best fit that minimize the distance
  dist_indexes.noise <- sort(dist_vals_noise, index.return=T)
  save_indexes.noise <- dist_indexes.noise$ix[1:ss-1]
  
  
  #choose 200 noisy samples to calculate the p values for all paramters
  
  noisy.sample=matrix(NA,nrow=length(noisy.observation),ncol=ss-1)
  
  
  
  for(kkk in 1:(ss-1)){
    noisy.sample[,kkk]=noisy.simulation.matrix[,save_indexes.noise[kkk]]
  }
  
  
  #the estimated theta + the best fit noisy data
  
  p1.vec=c(NA)
  p2.vec=c(NA)
  beta_0=c(NA)
  gamma_0=c(NA)
  
  for(jj in 1:(ss-1)){
    w0=noisy.sample[,jj]
    beta_0=saved_beta_vec[dist_indexes.noise$ix[jj]] 
    gamma_0=saved_gamma_vec [dist_indexes.noise$ix[jj]]
    remina.non.noisy=matrix(NA,nrow=length(noisy.observation),ncol=(ss-1)-1)
    ssss=1
    for(dd in 1:(ss-1)){
      
      if((dd != save_indexes.noise[jj])){
        
        remina.non.noisy[,ssss]=Accepted.simulation.matrix[,dd]
        ssss=ssss+1
      } 
      
    }
    
    
    ##########
    S_0=c(NA)
    
    for(iii in 1:(ss-1-1)){
      
      error.term_0=remina.non.noisy[,iii]-w0
      
      
      S_0[iii]=sum((error.term_0 /sd.e.term)^2)
      
    }
    
    
    accept.probability_0=c(NA)
    
    for(iii in 1:length(S_0)){
      accept.probability_0[iii]= dchisq(S_0[iii],length(noisy.observation) )*sqrt(S_0[iii])
    }
    
    c_0=min(accept.probability_0) 
    
    normlize.accept.probability_0=c(NA)
    
    for(iii in 1:length(S_0)){
      normlize.accept.probability_0[iii]= accept.probability_0[iii]/c_0
    }  
    
    Accepted.normlize.accept.probability_0=c(NA)
    saved_theta1_0=c(NA)
    saved_theta2_0=c(NA)
    
    gg=1
    one.coun=1
    two.coun=1
    
    Accepted.simulation.matrix_0=matrix(NA,nrow=length(noisy.observation),ncol=length(S_0))
    Accepted.probability.theta1=c(NA)
    Accepted.probability.theta2=c(NA)
    for(ddd in 1:length(S_0)){
      
      if(normlize.accept.probability_0[ddd]>=runif(1)){
        
        Accepted.normlize.accept.probability_0[gg]= normlize.accept.probability_0[ddd]
        
        if(saved_beta_vec[ddd]>=beta_0){
          saved_theta1_0[one.coun] <- saved_beta_vec[ddd]
          Accepted.probability.theta1[one.coun]= normlize.accept.probability_0[ddd]
          one.coun=one.coun+1
        }
        
        
        if(saved_gamma_vec[ddd]>=gamma_0){
          saved_theta2_0[two.coun] <- saved_gamma_vec[ddd]
          Accepted.probability.theta2[two.coun]= normlize.accept.probability_0[ddd]
          
          two.coun=two.coun+1
        }
        gg=gg+1
        
        
      }
      
    }
    
    
    
    
    #count p
    p1=sum(Accepted.probability.theta1)/sum(Accepted.normlize.accept.probability_0)
    
    p2=sum(Accepted.probability.theta2)/sum(Accepted.normlize.accept.probability_0)
    
    #p3=sum(Accepted.probability.theta3)/sum(Accepted.normlize.accept.probability_0)
    
    p1.vec[jj]=p1
    p2.vec[jj]=p2
    
    
    
  }
  
  output.df<- data.frame(p1.vec.beta=p1.vec, p2.vec.gamma=p2.vec, beta.samples=saved_beta_vec, gamma.samples=saved_gamma_vec)
  
  return(output.df) 
  
}

  
  

