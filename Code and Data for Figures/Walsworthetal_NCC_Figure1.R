#===============================================================================
#
# Plot effect of genetic variance, dispersal and protection level on ending
#   coral cover
#
#===============================================================================
setwd("C:/Users/Tim Walsworth/Dropbox/AWSCoralRuns")

iterout<-read.csv("DataforFigure1AB.csv",header=T)

nV<-length(unique(iterout$V))
nRS<-length(unique(iterout$ReserveSize))

library(gplots)

cols=rev(c("dodgerblue",NA,"darkorange",NA,"darkgreen"))

pdf("Walsworthetal_Fig1_Proofs.pdf",height=4,width=8)
layout(matrix(c(1,2,3),nrow=1,ncol=3,byrow=T),widths=c(1.1,1,.2))
par(cex.lab=1.3,xpd=T,mar=c(5,4,2,0)+.1)
plot(seq(0,6,length.out=8),rep(0,8),ylim=c(0,.7),xlim=c(0,6),xaxs="r",
     yaxs="i",type="n",ylab="Proportion Coral Cover",
     xlab=expression(paste("Genetic Variance (",V[i],")",sep="")),bty="n",xaxt="n",
     yaxt="n")
axis(side=1,at=seq(0.5,5.5,by=1),labels=c(unique(iterout$V)[order(unique(iterout$V))][-2]))
axis(side=2,at=c(0,.15,.3,.45,.6),las=2)
randonly<-iterout[iterout$Strategy==9,]
sub1<-randonly[randonly$D==0,]


for(i in 1:nRS)
{
  if(i ==2 |i == 4) next
  subv.all<-sub1[sub1$ReserveSize==unique(iterout$ReserveSize)[order(unique(iterout$ReserveSize))][i],]
  
  
  subv<-aggregate(subv.all$AvDensity,list("V"=subv.all$V),quantile,probs=c(.05,.25,.5,.75,.95))[-2,]
  print(subv)
  
  segments(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,2],seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,4],col=cols[i],lwd=4)
  segments(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,1],seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,5],col=cols[i],lwd=1)
  points(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,3],pch="_",col="black",cex=1.5)
 
}
legend("topleft","No Dispersal",bty="n")


par(mar=c(5,0.5,2,.2))

plot(seq(0,6,length.out=8),rep(0,8),ylim=c(0,.7),xlim=c(0,6),xaxs="i",
     yaxs="i",type="n",ylab="Coral Cover",
     xlab=expression(paste("Genetic Variance (",V[i],")",sep="")),
     bty="n",xaxt="n",
     yaxt="n")
axis(side=1,at=seq(0.5,5.5,by=1),labels=c(unique(iterout$V)[order(unique(iterout$V))][-2]))

randonly<-iterout[iterout$Strategy==9,]
sub1<-randonly[randonly$D==.01,]

for(i in 1:nRS)
{
  if(i ==2 |i == 4) next
  subv.all<-sub1[sub1$ReserveSize==unique(iterout$ReserveSize)[order(unique(iterout$ReserveSize))][i],]
  
  
  subv<-aggregate(subv.all$AvDensity,list("V"=subv.all$V),quantile,probs=c(.05,.25,.5,.75,.95))[-2,]
  print(subv)

  segments(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,2],seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,4],col=cols[i],lwd=4)
  segments(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,1],seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,5],col=cols[i],lwd=1)
  points(seq(0.5,5.5,by=1)+((-1)^i)*.1+(i-1)*.05,subv[,2][,3],pch="_",col="black",cex=1.5)
  
  
}

legend("topleft","High Dispersal",bty="n")
par(mar=c(0,0,0,1))
frame()
legend("center",lty=1,col=cols[c(5,3,1)],
       legend=unique(iterout$ReserveSize)[rev(order(unique(iterout$ReserveSize)))][c(1,3,5)],
       lwd=2,y.intersp=1.5, horiz=F,bty="n",title="Proportion\nProtected")


dev.off()

