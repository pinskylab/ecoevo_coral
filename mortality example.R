source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1)
scenarios<-read.csv("Norberg_StartingParameters_scenariolist.csv",header=T)

aa<-seq(20,30,by=.1)
bb<-mortality.fun(TC=aa,zi=25,rmax=1.5,w=1.5,growth="normal")
cc<-growth.fun(rmax=1.5,TC=aa,zi=25,w=1.5,growth="normal")
dd<-growth.fun(rmax=1.5,TC=aa,zi=25,w=4.5,growth="normal")
ee<-mortality.fun(TC=aa,zi=25,rmax=1.5,w=4.5,growth="normal")
plot(aa,cc,ylim=c(-1,1.5),type="l",lty=2,ylab="Rate",
     xlab="Temperature",main="Growth, mortality, and lambda for Topt=25")
points(aa,bb,lty=1,type="l")
points(aa,dd,lty=1,type="l",col="red")
points(aa,ee,lty=2,col="red",type="l")
lines(aa,(cc-bb),lty=1,lwd=2)
lines(aa,(dd-ee),lty=1,lwd=2,col="red")
abline(h=0,lty=3)
# legend("top",lty=c(2,1,1),lwd=c(1,1,2),col=c("black","black","black"),
#        legend=c("Competitive Growth","Competitive Mortality","Competitive Lambda"),bty="n")

legend("top",lty=c(2,1,1,2,1,1),lwd=c(1,1,2,1,1,2),col=c("black","black","black","red","red","red"),
       legend=c("Competitive Growth","Competitive Mortality","Competitive Lambda",
                "Tolerant Growth","Tolerant Mortality","Tolerant Lambda"),bty="n")
abline(v=c(25,26,27))
