# Make Figure 2 of Manuscript
# - Two panels
# - Top panel: Total coral cover on Y, proportion protected on X
#    - Different lines for the management strategies
#    - Intermediate dispersal and gen varaiance
# - Bottom panel: Relative coral cover boxplots
#    - Horizontal boxes
#    - Legend built into each line?
#====================================================================
load("Walsworthetal_NCC_Figure3Data.Rdata.rdata")
res_size<-0.2
jittervector<-seq(-.02,.02,length.out=8)
jittervector<-c(jittervector[1:7],NA,jittervector[8])
cols1<-c("lightblue","dodgerblue","dodgerblue4")
cols2<-colorRampPalette(c("burlywood","darkorange"))(2)
cols3<-colorRampPalette(c("green2","darkgreen"))(3)


cols<-c(cols1,cols2,cols3[c(1,2)],"black",cols3[3])
pdf("Walsworthetal_Fig3_Proofs.pdf",height=4,width=4,pointsize=4)
layout(matrix(c(1,1,3,
                1,1,3,
                2,2,3,
                2,2,4),nrow=4,ncol=3,byrow=T),widths=c(1,1,1.4))
par(xpd=T,mar=c(4,5,1,2))
boxsubdat<-iterout2[iterout2$Scenario==5 & iterout2$ReserveSize==res_size,]
my.boxplot.formula(boxsubdat$AvDensity~boxsubdat$Strategy,col=cols[-8],frame=F,border="black",
                   ylim=c(0,1),yaxt="n",horizontal=T,pch=16,xlab="Relative Coral Cover",par(fg="grey"),xaxs="i",yaxs="i",lwd=.5,cex=.7,outline=F)

for(i  in 1:(nstrat))
{
  if(i==8) next
  h<-i
  if(i==9) h<-8
  points(0.05,h,pch=16,cex=1.2,col=cols[i])
  lines(c(0.02,.08),c(h,h),col=cols[i],lwd=1)
  text(.08,h,stratlabs2[i],cex=1.5,pos=4,col="black")
  
}
text(0.035,8.7,"(a)",col="black",cex=1.8)

par(xpd=TRUE,mar=c(4,5,0,2),cex.lab=1.6)
plot(NA,NA,xlim=c(0,.5),ylim=c(0,.7),yaxs="i",xaxs="i",bty="l",
     ylab="Proportion Coral Cover",xlab="Proportion Protected",fg="grey")

for(i in 1:nstrat)
{
  if(i==8) next
  subdat<-iterout[iterout$Strategy==i & iterout$Scenario==5,]
  quants<-aggregate(subdat$AvDensity,list("Reserve"=subdat$ReserveSize),quantile,c(.05,.5,.95))
  if(i==9) abline(h=quants$x[3,2],col=cols[i],lty=2,xpd=F,lwd=.5)
  points(quants$Reserve+c(0,rep(jittervector[i],4)),quants$x[,2],lwd=1,col=c("black",rep(cols[i],4)),type="p",pch=16,cex=1.2)
  arrows(quants$Reserve+c(0,rep(jittervector[i],4)),quants$x[,1],quants$Reserve+c(0,rep(jittervector[i],4)),
           quants$x[,3],col=c("black",rep(cols[i],4)),length=.01,angle=90,code=3,lwd=1)

}
text(.025,.65,"(b)",col="black",cex=1.8)

plotorder<-rev(c(6,8,7,4,5,3,1,2))

stratlabs<-c("Hot Reefs","Cold Reefs","Hot and Cold\nReefs","High Cover","Low Cover","Evenly-spaced","Portfolio","Random")

cols<-colorRampPalette(rev(c("dodgerblue","darkred")))(100)
colbreaks<-seq(.01,1,by=.01)

par(mar=c(6,8,2,2),cex.main=1.5,cex.axis=1.5,cex.lab=1.5)

plot(NA,NA,ylim=c(0,8),xlim=c(0,3),yaxs="i",xaxs="i",ylab="",xlab="V",axes=F)
axis(side=2, at=seq(0.5,7.5,by=1),labels=stratlabs[plotorder],las=2,cex.axis=1.5)
axis(side=1,at=seq(0.5,2.5,by=1),labels=c(0,0.1,0.4))

ctr<-1
for(i in plotorder)
{
  
  for(j in 1:3)
  {
    
    xtable<-list(relcov1all,relcov2all,relcov3all)[[j]]
    max.yr<-apply(xtable,MARGIN=1,FUN=max)
    
    relval<-apply(xtable,MARGIN=2,FUN=function(x) x/max.yr)
    
    minrel<-apply(relval,MARGIN=1,min)
    
    relval1<-(relval[,i]-minrel)/(apply(relval,MARGIN=1,max)-minrel)
    colindex<-which.min(abs(colbreaks-mean(relval1,na.rm=T)))
   
    if(colbreaks[colindex]-mean(relval1,na.rm=T)<0)colindex<-colindex+1
    polygon(c(j-1,j,j,j-1),c(ctr-1,ctr-1,ctr,ctr),col=adjustcolor(cols[colindex],
                                                                  alpha.f=1),border=F)
    
  }
  ctr<-ctr+1
}

text(-.5,8.15,xpd=T,col="black","(c)",cex=1.8)

par(mar=c(6,8,2,2))
cols<-colorRampPalette(rev(c("dodgerblue","darkred")))(10)
plot(NA,NA,ylim=c(0,1),xlim=c(0,3),yaxs="i",xaxs="i",ylab="",xlab="",axes=F,bty="o",main="Relative Performance")
axis(side=1,at=c(0,3),labels=c("Low","High"))

for(i in 1:10){

    polygon(c((i-1)*.3,i*.3,i*.3,(i-1)*.3),c(0,0,1,1),
              col=adjustcolor(cols[i],alpha.f=1),border=F)
    
}

dev.off()

