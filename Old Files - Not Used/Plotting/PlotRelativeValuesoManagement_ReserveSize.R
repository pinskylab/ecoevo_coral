rm(list=ls())
graphics.off()
load("C:/Users/Tim Walsworth/Dropbox/AWSEcoEvoCoral/Base_ReserveSize_CoralSimulationsRuns010317.RData")
scenarios<-read.csv("Norberg_StartingParameters_iterationlist_reservesize.csv",header=T)
strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","portfolioGreedy","random") # Name different management strategies
nstrat<-length(strategies)
reservesize<-unique(scenarios$reserve)[1:5]
#which(scenarios$reserve==reservesize[3])
# zero20<-read.csv("ReserveSizeSimulations_ZeroTwentyProtection_new.csv",header=T)
# zero<-zero20[zero20[,3]==0,]
# twenty<-zero20[zero20[,3]==0.2,]
# for(i in 1:6)
# {
#   zero1<-zero
#   zero1[,2]<-i
#   if(i ==1) zeroall<-zero1
#   else zeroall<-rbind(zeroall,zero1)
# }

library(RColorBrewer)
str(tsoutRuns)

niter<-iterations
modname<-"Base_Heatmaps"
monitortime<-100
monitor.yrs<-c(20,50,100,500)
size<-60
niter<-50
# 1. For each scenario and each iteration, calculate average coral cover across reef
# 2. Rank management strategies
# 3. Calculate relative values of different strategies
# 4. Plot average relative values (barplot with error bars?)
rankcorals<-array(NA,dim=c(niter,nstrat,9,length(reservesize)))
relcorals<-array(NA,dim=c(niter,nstrat,9,length(reservesize)))
maxavcorals<-array(NA,dim=c(niter,9,length(reservesize)))
logratcorals<-array(NA,dim=c(niter,nstrat,9,length(reservesize)))
wtvartraitsp1<-logratcorals
wtvartraitsp2<-logratcorals
extinctrate1<-logratcorals
extinctrate2<-extinctrate1
extinctrateboth<-extinctrate1
avcorals<-relcorals

wtd.relcorals<-matrix(NA,nrow=6*niter,ncol=nstrat)

for(i in 1:length(unique(scenarios$Scenario)))
{

  for(j in 1:50)
  {

    for(z in 1:length(reservesize))
    {
    print(i)
      print(z)
      # i=1
      # j=1
      # z=1
    if((i+9*(z-1))<=16) subdat<-tsoutRuns[[(i+9*(z-1))]][1:(7*size),which(monitor.yrs==monitortime)+1,j,]
    if((i+9*(z-1))>16 & (i+9*(z-1)) <=32) subdat<-tsoutRunssecond[[(i+9*(z-1))-16]][1:(7*size),which(monitor.yrs==monitortime)+1,j,]
    if((i+9*(z-1))>32) subdat<-tsoutRunsThird[[(i+9*(z-1))-32]][1:(7*size),which(monitor.yrs==monitortime)+1,j,]
      # subdat<-tsoutRuns[[(i+9*(z-1))]][1:(7*size),which(monitor.yrs==monitortime)+1,j,]
      
      sp1<-subdat[1:size,]
    #print(sp1)
    sp2<-subdat[(size+1):(2*size),]
    allcorals<-subdat[1:(2*size),]
    trait1<-subdat[(3*size+1):(4*size),]
    trait2<-subdat[(4*size+1):(5*size),]
    alltrait<-subdat[(3*size+1):(5*size),]
    temps<-subdat[(6*size+1):(7*size),]
    corals<-sp1+sp2
    #print(corals)
    avcorals[j,,i,z]<-colSums(corals)/size
    maxavcorals[j,i,z]<-max(avcorals[j,,i,z])
    for(k in 1:length(strategies))
    {
      if((i+9*(z-1))<=16){
        extinctrate1[j,k,i,z]<-tsoutRuns[[(i+9*(z-1))]][1,4*length(monitor.yrs)+2,j,k]
        extinctrate2[j,k,i,z]<-tsoutRuns[[(i+9*(z-1))]][1,4*length(monitor.yrs)+3,j,k]
        extinctrateboth[j,k,i,z]<-tsoutRuns[[(i+9*(z-1))]][1,4*length(monitor.yrs)+4,j,k]
      }
      if((i+9*(z-1))>16 & (i+9*(z-1)) <=32){
        extinctrate1[j,k,i,z]<-tsoutRunssecond[[(i+9*(z-1))-16]][1,4*length(monitor.yrs)+2,j,k]
        extinctrate2[j,k,i,z]<-tsoutRunssecond[[(i+9*(z-1))-16]][1,4*length(monitor.yrs)+3,j,k]
        extinctrateboth[j,k,i,z]<-tsoutRunssecond[[(i+9*(z-1))-16]][1,4*length(monitor.yrs)+4,j,k]
      }
      if((i+9*(z-1))>32){
        extinctrate1[j,k,i,z]<-tsoutRunsThird[[(i+9*(z-1))-32]][1,4*length(monitor.yrs)+2,j,k]
        extinctrate2[j,k,i,z]<-tsoutRunsThird[[(i+9*(z-1))-32]][1,4*length(monitor.yrs)+3,j,k]
        extinctrateboth[j,k,i,z]<-tsoutRunsThird[[(i+9*(z-1))-32]][1,4*length(monitor.yrs)+4,j,k]
      }
      # extinctrate1[j,k,i,z]<-tsoutRuns[[(i+9*(z-1))]][1,4+which(monitor.yrs==monitortime)+1,j,k]
      # extinctrate2[j,k,i,z]<-tsoutRuns[[(i+9*(z-1))]][1,8+which(monitor.yrs==monitortime)+1,j,k]
      # # wtvartraitsp1[j,k,i,z]<-wtd.var(x=trait1[,k],weights=sp1[,k])  # weighted variance of trait sp1
      # wtvartraitsp2[j,k,i,z]<-wtd.var(x=trait2[,k],weights=sp2[,k])  # weighted variance of trait sp2
      wtvartraitsp1[j,k,i,z]<-sum((alltrait[,k]*allcorals[,k])/sum(allcorals[,k]))-sum((rep(temps[,k],2)*allcorals[,k])/sum(allcorals[,k])) # weighted difference between coral traits and environment
      wtvartraitsp2[j,k,i,z]<-sum((trait2[,k]*sp2[,k])/sum(sp2[,k]))-sum((temps[,k]*sp2[,k])/sum(sp2[,k]))
    } 
    }
    mac<-max(maxavcorals[j,i,])
    for(z in 1:length(reservesize))
    {
      relcorals[j,,i,z]<-avcorals[j,,i,z]/mac
      rankcorals[j,,i,z]<-rank(avcorals[j,,i,z])
      logratcorals[j,,i,z]<-log(colMeans(sp1)/colMeans(sp2))
      
      
    }
    
  }
  
}

iterout<-matrix(NA,nrow=length(reservesize)*length(unique(scenarios$Scenario))*niter*nstrat,ncol=13)
colnames(iterout)<-c("Scenario","Strategy","ReserveSize","Iteration","AvDensity",
                     "LogRatCorals","WVTraitSp1","WVTraitSp2","Ext1","Ext2","ExtBoth","V","D")

ctr<-1
for(j in 1:niter)
{
  for(i in 1:length(unique(scenarios$Scenario)))
  {
    for(k in 1:length(strategies))
    {
      for(z in 1:length(reservesize))
      {
        iterout[ctr,1]<-as.numeric(i)
        iterout[ctr,2]<-k
        iterout[ctr,3]<-as.numeric(reservesize[z])
        iterout[ctr,4]<-as.numeric(j)
        iterout[ctr,5]<-as.numeric(avcorals[j,k,i,z])
        iterout[ctr,6]<-as.numeric(logratcorals[j,k,i,z])
        iterout[ctr,7]<-as.numeric(wtvartraitsp1[j,k,i,z])
        iterout[ctr,8]<-as.numeric(wtvartraitsp2[j,k,i,z])
        iterout[ctr,9]<-as.numeric(extinctrate1[j,k,i,z])
        iterout[ctr,10]<-as.numeric(extinctrate2[j,k,i,z])
        iterout[ctr,11]<-as.numeric(extinctrateboth[j,k,i,z])
        iterout[ctr,12]<-scenarios[1:9,2][which(scenarios[,1]==iterout[ctr,1])[1]]
        iterout[ctr,13]<-scenarios[1:9,3][which(scenarios[,1]==iterout[ctr,1])[1]]
        ctr<-ctr+1
      }
    }
    
  }
  
}


write.csv(iterout,file=paste(modname,monitortime,"ReserveSizeSimulations_AverageDensities_newStrats_",Sys.Date(),".csv",sep=""),row.names = F)

iterout2<-iterout
reservesize<-unique(iterout[,3])[order(unique(iterout[,3]))]
for(i in 1:length(unique(scenarios$Scenario)))
{
  for(k in reservesize)
  {
    
 
  dat<-iterout[iterout[,1]==i & iterout[,3]==k,]
  
  for(j in 1:niter)
  {
    whichscen<-which(iterout[,1]==i&iterout[,4]==j & iterout[,3]==k,arr.ind=T)
    subdat<-dat[dat[,4]==j,]
    maxdens<-max(subdat[,5])
    iterout2[whichscen,5]<-subdat[,5]/maxdens
  }
  }
}



write.csv(iterout2,file=paste(modname,monitortime,"ReserveSizeSimulations_RelativeDensities_newStrats_",Sys.Date(),".csv",sep=""),row.names = F)
#iterout<-read.csv("ReserveSizeSimulations_AverageDensities_newStrats.csv",header=T)
#iterout2<-read.csv("ReserveSizeSimulations_RelativeDensities_Test.csv")

# subdat<-iterout[iterout[,1]==7,]
# 
# a<-aggregate(subdat[,5],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.05,.10,.25,.50,.75,.90,.95))
# jittervector<-c(-.005,.005,-.005,.005,-.005,.005)
# plot(seq(-.050,.55,length.out=10),seq(0,1,length.out=10),type="n",ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
# for(i in 1:nstrat)
# {
#   subdat2<-a[a$Strat==i,]
#   #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[i],alpha.f=.5))
#   lines(subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,4],col=cols[i],lwd=2)
#   arrows(subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,4],subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,7],
#          col=c("black",rep(cols[i],6)),angle=90,length=.001)
#   arrows(subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,4],subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,1],
#          col=c("black",rep(cols[i],6)),angle=90,length=.001)
#   points(subdat2$Size+jittervector[i]*c(0,rep(1,6)),subdat2$x[,4],col=c("black",rep(cols[i],6)),pch=16,cex=2)
# }
# legend("topleft",bty="n",col=cols,legend=strategies,cex=3,horiz=F,pch=15)

#cols<-brewer.pal(7,"Dark2")
lreserve<-length(reservesize)

cols<-brewer.pal(9,"Paired")
pdf1name<-paste(modname,"_RelativeDensitiesManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
#windows()
pdf(pdf1name,height=12,width=18)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,7,5,5))
jittervector<-c(-.005,.015,-.01,.01,-.015,.005,-.02,.02,0)
for(i in 1:9)
{
  subdat<-iterout2[iterout2[,1]==i,]
  # print(subdat)
  a<-aggregate(subdat[,5],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
  par(mar=c(.1,.1,.1,.1),bg="white")
  plot(seq(-.01,.8,length.out=10),seq(0,1.05,length.out=10),type="n",ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
  
  if(i==7){
    mtext(side=1,"Proportion Protected",line=2.5)
    mtext(side=2,"Relative Cover",line=2.5)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,2]
    #        ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,8]
    #        ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],lreserve-1)),pch=16,cex=2)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    #lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    arrows(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,2]
           ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    arrows(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,8]
           ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    points(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],col=c("black",rep(cols[j],lreserve-1)),pch=16,cex=2)
  }
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=2,horiz=F,ncol=3,pch=15)
dev.off()


#cols<-brewer.pal(7,"Dark2")
pdf2name<-paste(modname,"_AverageDensitiesManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
#windows()
pdf(pdf2name,height=12,width=18)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,7,5,5))
jittervector<-c(-.005,.015,-.01,.01,-.015,.005,-.02,.02,0)*1.5
for(i in 1:9)
{
  subdat<-iterout[iterout[,1]==i,]
  
  a<-aggregate(subdat[,5],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
  par(mar=c(.1,.1,.1,.1),bg="white")
  plot(seq(-.01,.8,length.out=10),seq(0,1.05,length.out=10),type="n",ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
  
  if(i==7){
    mtext(side=1,"Proportion Protected",line=2.5)
    mtext(side=2,"Average Cover",line=2.5)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,2]
    #        ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,8]
    #        ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],lreserve-1)),pch=16,cex=2)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    #lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    arrows(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,2]
           ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    arrows(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,8]
           ,col=c("black",rep(cols[j],lreserve-1)),angle=90,length=.001)
    points(subdat2$Size+jittervector[j]*c(0,rep(1,lreserve-1)),subdat2$x[,5],col=c("black",rep(cols[j],lreserve-1)),pch=16,cex=2)
  }
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=F,ncol=3,pch=15)
dev.off()



example.size<-0.3
example.dat<-iterout[iterout[,3]==example.size,]
example.dat<-cbind(example.dat,example.dat[,5])

for(i in 1:9)
{
  for(j in 1:50)
  {
    subdat<-example.dat[example.dat[,1]==i&example.dat[,4]==j,]
    example.dat[which(example.dat[,1]==i&example.dat[,4]==j,arr.ind = T),6]<-example.dat[which(example.dat[,1]==i&example.dat[,4]==j,arr.ind = T),5]/
      max(subdat[,5])
    

    
  }
}



cols<-brewer.pal(9,"Paired")
pdf1name<-paste(modname,"_RelativeDensitiesManagementStrategiesBoxplot_",monitortime,"_",example.size,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
#windows()
pdf(pdf1name,height=12,width=18)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,5,5,5))
for(i in 1:9)
{
  subdat<-example.dat[example.dat[,1]==i,]
  par(mar=c(.1,.1,.1,.1),bg="white")
  a<-boxplot(subdat[,6]~subdat[,2],xlim=c(.5,9.5),col=cols,yaxs="i",ylim=c(0,1),frame=F,
          yaxt="n",horizontal=T)
  print(a$stats)
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
legend("center",bty="n",col=cols,legend=strategies,cex=2,horiz=F,ncol=3,pch=15)
dev.off()

# 
# 
# #cols<-brewer.pal(7,"Dark2")
# pdf3name<-paste(modname,"_LogRatCoralSppManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
# #windows()
# pdf(pdf3name,height=12,width=18)
# layout(matrix(c(1,0,2,0,3,
#                 0,0,0,0,0,
#                 4,0,5,0,6,
#                 0,0,0,0,0,
#                 7,0,8,0,9,
#                 10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
# par(oma=c(0,7,5,5))
# jittervector<-c(-.005,.015,-.01,.01,-.015,.005)
# for(i in 1:9)
# {
#   subdat<-iterout[iterout[,1]==i,]
#   
#   a<-aggregate(subdat[,6],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
#   par(mar=c(.1,.1,.1,.1),bg="white")
#   plot(seq(-.01,.55,length.out=10),seq(-15,15,length.out=10),type="n",
#        ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
#   
#   if(i==7){
#     mtext(side=1,"Proportion Protected",line=2.5)
#     mtext(side=2,"log(Comp/ST)",line=2.5)
#   }
#   for(j in 1:nstrat)
#   {
#     subdat2<-a[a$Strat==j,]
#     #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
#     # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
#     # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
#     #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
#     #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
#     
#     lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
#     arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,1]
#            ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,9]
#            ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
#   }
# }
# mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
# mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
# mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
# mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
# mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
# mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
# frame()
# legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
# dev.off()
# 
# #cols<-brewer.pal(7,"Dark2")
# pdf4name<-paste(modname,"_TraitDiffCoralsManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
# #windows()
# pdf(pdf4name,height=12,width=18)
# layout(matrix(c(1,0,2,0,3,
#                 0,0,0,0,0,
#                 4,0,5,0,6,
#                 0,0,0,0,0,
#                 7,0,8,0,9,
#                 10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
# par(oma=c(0,7,5,5))
# jittervector<-c(-.005,.015,-.01,.01,-.015,.005)
# for(i in 1:9)
# {
#   subdat<-iterout[iterout[,1]==i,]
#   
#   a<-aggregate(subdat[,7],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
#   par(mar=c(.1,.1,.1,.1),bg="white")
#   plot(seq(-.01,.55,length.out=10),seq(-8,5,length.out=10),type="n",
#        ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
#   
#   if(i==7){
#     mtext(side=1,"Proportion Protected",line=2.5)
#     mtext(side=2,expression(paste(Trait[w]-Temp[w])),line=2.5)
#   }
#   for(j in 1:nstrat)
#   {
#     subdat2<-a[a$Strat==j,]
#     #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
#     # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
#     # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
#     #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
#     #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
#     
#     lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
#     arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,1]
#            ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,9]
#            ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
#     points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
#   }
#   abline(h=1.2,lty=2,lwd=1)
# }
# mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
# mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
# mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
# mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
# mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
# mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
# frame()
# legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)
# dev.off()
# 
extinctval<-.01
whichspp<-1
pdf1name<-paste(modname,"_ExtinctionsThroughTimeManagementStrategies_",example.size,"_",extinctval,"_",ifelse(whichspp==1,"Comp",ifelse(whichspp==2,"Tol","Both")),"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
pdf(pdf1name,height=4,width=6,pointsize=8)
# layout(matrix(c(1,2,
#                 3,4,
#                 5,6,
#                 7,8),nrow=2,ncol=4,byrow=T))
layout(matrix(c(1,2,3,4,5,6,
                7,8,9,10,11,12,
                13,14,15,16,17,18),nrow=3,ncol=6,byrow=T))
par(mar=c(5,4,0,0),oma=c(2,4,4,2))
for(k in 1:9)
{
  subdat<-iterout[iterout[,1]==k &iterout[,3]==.30,]
  subdat[is.na(subdat[,9]),9]<-501
  subdat[is.na(subdat[,10]),10]<-501
  subdat[is.na(subdat[,11]),11]<-501
  outdat<-matrix(NA,nrow=500,ncol=9)
  for(i in 1:9)
  {
    for(j in 1:500)
    {
      subdat2<-subdat[subdat[,2]==i,]
      outdat[j,i]<-sum(subdat2[,8+whichspp]<=j)
    }
  }
  plot(seq(1,500,by=1),outdat[,1],col=cols[1],type="l",lwd=1,ylim=c(0,50),ylab="Extinctions",xlab="Year",bty="l",las=2)
  for(i in 2:9)
  {
    lines(seq(1,500,by=1),outdat[,i],col=cols[i],lwd=1)
  }
  barplot(outdat[j,],col=cols,ylim=c(0,50),border=F,ylab="Extinctions/50 iterations",names=strategies,las=2)
  
  
}

mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=2)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=2)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=2)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)


# mtext(side=2,line=2,at=.85,"D= 0.001",outer=T,cex=1.5)
# mtext(side=2,line=2,at=.35,"D= 0.01",outer=T,cex=1.5)
# mtext(side=3,line=2,at=.25,"V=0.1",outer=T,cex=1.5)
# mtext(side=3,line=2,at=.75,"V=0.4",outer=T,cex=1.5)
dev.off()


pdf1name<-paste(modname,"_PersistenceManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
#windows()
pdf(pdf1name,height=12,width=18)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,7,7,5),cex.axis=1.5)
for(i in 1:9)
{
  par(mar=c(.1,.1,.1,.1),bg="white")
  plot(seq(0,.8,length.out=10),seq(0,1.01,length.out=10),type="n",bty="l",
       ylab=ifelse(i==7,"Persistence",""),xlab=ifelse(i==7,"Proportion MPA",""),yaxs="i")
  if(i==7){
    mtext(side=1,"Proportion Protected",line=2.5)
    mtext(side=2,"Persistence",line=2.5)
  }
  for(z in 1:9)
  {
    dat.out<-matrix(NA,nrow=lreserve,ncol=2)
    for(j in 1:lreserve)
    { 
      
      #i=1;j=1;z=1
      subdat<-iterout[iterout[,1]==i &iterout[,3]==reservesize[j]&iterout[,2]==z,]
      dat.out[j,1]<-reservesize[j]
      dat.out[j,2]<-1-(sum(!is.na(subdat[,11]))/50)
      
    
    }
    lines(dat.out[,1],dat.out[,2],col=cols[z],lwd=2)
  }
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=2,horiz=F,ncol=3,lty=1,lwd=2)
dev.off()



#cols<-brewer.pal(7,"Dark2")
pdf5name<-paste(modname,"_ExtinctRateCoralsManagementStrategies_ReserveSize_lines",monitortime,"_",format(Sys.Date(),format="%m%d%y"),".pdf",sep="")
#windows()
pdf(pdf5name,height=12,width=18,onefile=T)
layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,7,5,5))
jittervector<-c(-.005,.015,-.01,.01,-.015,.005)
for(i in 1:9)
{
  subdat<-iterout[iterout[,1]==i,]
  
  a<-aggregate(subdat[,11],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
  par(mar=c(.1,.1,.1,.1),bg="white")
  plot(seq(-.01,.55,length.out=10),seq(0,1.05,length.out=10),type="n",
       ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
  
  if(i==7){
    mtext(side=1,"Proportion Protected",line=2.5)
    mtext(side=2,expression(paste(Trait[w]-Temp[w])),line=2.5)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,3]
           ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,7]
           ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
  }
  abline(h=1.2,lty=2,lwd=1)
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)

layout(matrix(c(1,0,2,0,3,
                0,0,0,0,0,
                4,0,5,0,6,
                0,0,0,0,0,
                7,0,8,0,9,
                10,10,10,10,10),nrow=6,ncol=5,byrow=T),heights=c(1,.3,1,.3,1,.5),widths=c(1,.2,1,.2,1))
par(oma=c(0,7,5,5))
jittervector<-c(-.005,.015,-.01,.01,-.015,.005)
for(i in 1:9)
{
  subdat<-iterout[iterout[,1]==i,]
  
  a<-aggregate(subdat[,10],list("Strat"=subdat[,2],"Size"=subdat[,3]),quantile,probs=c(.01,.05,.10,.25,.50,.75,.90,.95,.99))
  par(mar=c(.1,.1,.1,.1),bg="white")
  plot(seq(-.01,.55,length.out=10),seq(0,1.05,length.out=10),type="n",
       ylab="Relative Hard Coral Density",xlab="Percent Protected",yaxs="i",xaxs="i")
  
  if(i==7){
    mtext(side=1,"Proportion Protected",line=2.5)
    mtext(side=2,expression(paste(Trait[w]-Temp[w])),line=2.5)
  }
  for(j in 1:nstrat)
  {
    subdat2<-a[a$Strat==j,]
    #polygon(c(subdat2$Size,rev(subdat2$Size)),c(subdat2$x[,2],rev(subdat2$x[,6])),col=adjustcolor(cols[j],alpha.f=.5))
    # lines(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=cols[j],lwd=2)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,1]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # arrows(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,9]
    #        ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    # points(subdat2$Size+jittervector[j]*c(0,rep(1,6)),subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
    
    lines(subdat2$Size+jittervector[j],subdat2$x[,5],col=cols[j],lwd=2)
    arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,3]
           ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    arrows(subdat2$Size+jittervector[j],subdat2$x[,5],subdat2$Size+jittervector[j],subdat2$x[,3]
           ,col=c("black",rep(cols[j],6)),angle=90,length=.001)
    points(subdat2$Size+jittervector[j],subdat2$x[,5],col=c("black",rep(cols[j],6)),pch=16,cex=2)
  }
  abline(h=1.2,lty=2,lwd=1)
}
mtext(side=2,at=.25,outer=T,"D = .01",cex=2,line=4)
mtext(side=2,at=.58,outer=T,"D = .001",cex=2,line=4)
mtext(side=2,at=.91,outer=T,"D = 0",cex=2,line=4)
mtext(side=3,at=.15,outer=T,"V = 0",cex=2,line=2)
mtext(side=3,at=.5,outer=T,"V = .1",cex=2,line=2)
mtext(side=3,at=.85,outer=T,"V = .4",cex=2,line=2)
frame()
legend("center",bty="n",col=cols,legend=strategies,cex=3,horiz=T,pch=15)

dev.off()


