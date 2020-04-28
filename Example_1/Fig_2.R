#plot Fig1
source("Initial_parameters_final.R", echo = TRUE)

plot(true_time,out[,3],typ="l",col=2,ylab="Proportion infected",xlab="Time(week)")
points(true_time[min_time:max_time], data.rep,pch=16,col=4)
legend("topright", legend = c("Observations", "True infection curve"), col = c("blue", "red"), lty=c(NA,1), pch=c(16,NA))
