aa<-seq(20,30,by=.1)
bb<-mortality.fun(TC=aa,zi=25,rmax=1.5,w=1,growth="normal")
cc<-growth.fun(rmax=1.5,TC=aa,zi=25,w=1,growth="normal")
dd<-growth.fun(rmax=1.5,TC=aa,zi=25,w=3,growth="normal")
ee<-mortality.fun(TC=aa,zi=25,rmax=1.5,w=3,growth="normal")
plot(aa,bb,ylim=c(-1,1.5),type="l",lty=2,ylab="Rate",
     xlab="Temperature",main="Growth, mortality, and lambda for Topt=25")
points(aa,cc,lty=1,type="l")
points(aa,dd,lty=1,type="l",col="red")
points(aa,ee,lty=2,col="red",type="l")
lines(aa,(cc-bb),lty=1,lwd=2)
lines(aa,(dd-ee),lty=1,lwd=2,col="red")
abline(h=0,lty=3)
legend("top",lty=c(1,2,1,1,2,1),lwd=c(1,1,2,1,1,2),col=c("black","black","black","red","red","red"),
       legend=c("Competitive Growth","Competitive Mortality","Competitive Lambda",
                "Tolerant Growth","Tolerant Mortality","Tolerant Lambda"),bty="n")
