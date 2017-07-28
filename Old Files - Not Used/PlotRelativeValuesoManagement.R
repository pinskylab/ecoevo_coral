#load("C:/Users/larogers-admin/Dropbox/AWSEcoEvoCoral/CoralSimulationsRuns102816.RData")
library(RColorBrewer)
str(tsoutRuns)

# 1. For each scenario and each iteration, calculate average coral cover across reef
# 2. Rank management strategies
# 3. Calculate relative values of different strategies
# 4. Plot average relative values (barplot with error bars?)
rankcorals<-array(NA,dim=c(50,7,9))
relcorals<-array(NA,dim=c(50,7,9))
maxavcorals<-matrix(NA,nrow=50,ncol=9)
logratcorals<-array(NA,dim=c(50,7,9))
for(i in 1:9)
{
  for(j in 1:50)
  {
    subdat<-tsoutRuns[[i]][1:120,1,j,]
    sp1<-subdat[1:60,]
    sp2<-subdat[61:120,]
    corals<-sp1+sp2
    avcorals<-colSums(corals)/60
    maxavcorals[j,i]<-max(avcorals)
    relcorals[j,,i]<-avcorals/max(avcorals)
    rankcorals[j,,i]<-rank(avcorals)
    logratcorals[j,,i]<-log(colMeans(sp1)/colMeans(sp2))
  }
  
}


cols<-brewer.pal(7,"Dark2")

#windows()
pdf("RelativeDensitiesManagementStrategiesBoxplot_text_hightempgradient100yrs.pdf",height=12,width=15)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,5,5,5))
for(i in 1:9)
{
  par(mar=c(.1,.1,.1,.1),bg="white")
    boxplot(relcorals[,rev(1:7),i],xlim=c(.5,7.5),col=rev(cols),yaxs="i",ylim=c(0,1),frame=F,
            yaxt="n",horizontal=T)
  text(.3,2,labels=paste(signif(mean(maxavcorals[,i]),3)," +- ",signif(sd(maxavcorals[,i]),3)),
       cex=1.5)
    #hist(relcorals[,j,i],xlim=c(0,1),breaks=c(seq(0,1,by=.0125)))
#   barplot(colMeans(relcorals[,,i]))
#   print(colMeans(rankcorals[,,i]))
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=2)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=2)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=2)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
dev.off()

#windows()
pdf("RelativeRankingManagementStrategies.pdf",height=12,width=15)
par(oma=c(0,4,4,0))
layout(matrix(c(1,0,8,0,15,
                2,0,9,0,16,
                3,0,10,0,17,
                4,0,11,0,18,
                5,0,12,0,19,
                6,0,13,0,20,
                7,0,14,0,21,
                0,0,0,0,0,
                22,0,29,0,36,
                23,0,30,0,37,
                24,0,31,0,38,
                25,0,32,0,39,
                26,0,33,0,40,
                27,0,34,0,41,
                28,0,35,0,42,
                0,0,0,0,0,
                43,0,50,0,57,
                44,0,51,0,58,
                45,0,52,0,59,
                46,0,53,0,60,
                47,0,54,0,61,
                48,0,55,0,62,
                49,0,56,0,63,
                64,64,64,64,64),nrow=24,ncol=5,byrow=T),widths=c(1,.1,1,.1,1),heights=c(rep(1,7),3,rep(1,7),3,rep(1,7),4))
for(i in 1:9)
{

  for(j in 1:7)
  {
    par(mar=c(.1,.1,.1,.1),bg="white")
    hist(rankcorals[,j,i],xlim=c(.5,7.5),breaks=c(seq(0.5,7.5,by=1)),xaxt="n",yaxt="n",
         col=cols[j],xlab="",main="",border=adjustcolor("black",alpha.f=0))
    #hist(relcorals[,j,i],xlim=c(0,1),breaks=c(seq(0,1,by=.0125)))
  }
  axis(side=1)
  #   barplot(colMeans(relcorals[,,i]))
  #   print(colMeans(rankcorals[,,i]))
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=2)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=2)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=2)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
dev.off()


pdf("RelativeDensitiesManagementStrategies.pdf",height=12,width=15)
par(oma=c(0,4,4,0))
layout(matrix(c(1,0,8,0,15,
                2,0,9,0,16,
                3,0,10,0,17,
                4,0,11,0,18,
                5,0,12,0,19,
                6,0,13,0,20,
                7,0,14,0,21,
                0,0,0,0,0,
                22,0,29,0,36,
                23,0,30,0,37,
                24,0,31,0,38,
                25,0,32,0,39,
                26,0,33,0,40,
                27,0,34,0,41,
                28,0,35,0,42,
                0,0,0,0,0,
                43,0,50,0,57,
                44,0,51,0,58,
                45,0,52,0,59,
                46,0,53,0,60,
                47,0,54,0,61,
                48,0,55,0,62,
                49,0,56,0,63,
                64,64,64,64,64),nrow=24,ncol=5,byrow=T),widths=c(1,.1,1,.1,1),heights=c(rep(1,7),3,rep(1,7),3,rep(1,7),4))
for(i in 1:9)
{
  
  for(j in 1:7)
  {
    par(mar=c(.1,.1,.1,.1),bg="white")
#     hist(rankcorals[,j,i],xlim=c(.5,7.5),breaks=c(seq(0.5,7.5,by=1)),xaxt="n",yaxt="n",
#          col=cols[j],xlab="",main="",border=adjustcolor("black",alpha.f=0))
    hist(relcorals[,j,i],xlim=c(0,1),breaks=c(seq(0,1,by=.0125)),xaxt="n",yaxt="n",
         col=cols[j],xlab="",main="",border=adjustcolor("black",alpha.f=0))
  }
  axis(side=1)
  #   barplot(colMeans(relcorals[,,i]))
  #   print(colMeans(rankcorals[,,i]))
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=2)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=2)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=2)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
dev.off()


# times<-seq(1,2000)
# temps<-c(rep(27,1500),rep(0,500))
# annual.temp.change<-.02
# maxtemp<-30
# for(i in 1501:2000)
# {
#   temps[i]<-annual.temp.change*temps[i-1]*(1-(temps[i-1]/maxtemp))+temps[i-1]  
#   
# }
# 
# temps2<-temps
# for(i in 1:2000)
# {
#   if(i<=1500)  temps2[i]<-temps[i]+anoms.burn[1,i,1]
#   if(i>1500) temps2[i]<-temps2[i]+anoms.runs[1,i-1500,1]
# }
# plot(times,temps,ylim=c(22,35),type="l",lwd=4,xlab="Year",ylab="Temperature",bty="l")
# lines(times,temps2,lty=2,lwd=2,col="blue")

#windows()
pdf("CoralSpLogRatsManagementStrategies.pdf",height=12,width=15)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,5,5,5))
for(i in 1:9)
{
  par(mar=c(.1,.1,.1,.1),bg="white")
  boxplot(logratcorals[,1:7,i],xlim=c(.5,7.5),col=cols,yaxs="i",yaxt="n",ylim=c(-15,15),xaxt="n",
          frame=F,horizontal=F)
  abline(h=0,lty=2)
  axis(side=2, at=c(-10,0,10),labels=c("ST","Even","C"))

}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=2)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=2)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=2)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
dev.off()



# Single panel as example
windows()

layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(4,1))

  par(mar=c(5,2,1,2),bg="white")
  boxplot(relcorals[,c(rev(1:6),7),8],xlim=c(.5,7.5),col=cols[c(rev(1:6),7)],yaxs="i",ylim=c(0,1),frame=F,
          yaxt="n",horizontal=T,xlab="Relative Coral Density",cex.lab=1.5,pch=16,outcol=adjustcolor("grey50",alpha.f=.5))
  text(.2,1,labels=paste(signif(mean(maxavcorals[,8]),3)," +- ",signif(sd(maxavcorals[,8]),3)),
       cex=1.5)

par(mar=c(1,0,1,0))
frame()
legend("top",bty="n",col=cols[c(7,1:3)],legend=strategies[c(7,1:3)],cex=1.5,horiz=T,pch=15,pt.cex=3)
legend("bottom",bty="n",col=cols[4:6],legend=strategies[4:6],cex=1.5,horiz=T,pch=15,pt.cex=3)


