# Make Figure 2 of Manuscript
# - Two panels
# - Top panel: Total coral cover on Y, proportion protected on X
#    - Different lines for the management strategies
#    - Intermediate dispersal and gen varaiance
# - Bottom panel: Relative coral cover boxplots
#    - Horizontal boxes
#    - Legend built into each line?
#====================================================================
scenarios<-read.csv("Norberg_StartingParameters_iterationlist_reservesize_SA.csv",header=T)
strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","portfolioGreedy","random") # Name different management strategies
nstrat<-length(strategies)
stratlabs<-c("Hot","Cold","Hot and Cold","High Cover","Low Cover","Gradient","Portfolio","Portfolio 2","Random")

iterout<-read.csv("Base_HeatmapsReserveSizeSimulations_AverageDensities_newStrats_2017-05-16.csv",header=T)
iterout2<-read.csv("Base_HeatmapsReserveSizeSimulations_RelativeDensities_newStrats_2017-05-16.csv",header=T)
iterout2<-iterout2[iterout2$Strategy!=8,]

res_size<-0.2
jittervector<-c(-.005,.015,-.01,.01,-.015,.005,-.02,.02,0)
cols<-brewer.pal(9,"Paired")

pdf("Walsworthetal_Fig2_CompareStrategies.pdf",height=4,width=3,pointsize=4)
layout(matrix(c(1,2),nrow=2,ncol=1))
par(xpd=TRUE,mar=c(4,4,1,2),cex.lab=1.2)
plot(NA,NA,xlim=c(0,.5),ylim=c(0,1),yaxs="i",xaxs="i",bty="l",
     ylab="Total Coral Cover",xlab="Proportion Protected")

for(i in 1:nstrat)
{
  if(i==8) next
  subdat<-iterout[iterout$Strategy==i & iterout$Scenario==5,]
  quants<-aggregate(subdat$AvDensity,list("Reserve"=subdat$ReserveSize),quantile,c(.05,.5,.95))
  points(quants$Reserve+c(0,rep(jittervector[i],4)),quants$x[,2],lwd=1,col=cols[i],type="o",pch=16,cex=1.2)
  segments(quants$Reserve+c(0,rep(jittervector[i],4)),quants$x[,1],quants$Reserve+c(0,rep(jittervector[i],4)),
           quants$x[,3],col=cols[i])
}
text(.015,.95,"(a)")
par(xpd=T,mar=c(4,4,0,2))
boxsubdat<-iterout2[iterout2$Scenario==5 & iterout2$ReserveSize==res_size,]
boxplot(boxsubdat$AvDensity~boxsubdat$Strategy,col=cols[-8],frame=F,
        ylim=c(0,1),yaxt="n",horizontal=T,pch=16,xlab="Relative Coral Cover",xaxs="i",yaxs="i",lwd=.5,cex=.8)

#polygon(c(0,.4,.4,0),c(.3,.3,8.7,8.7),col="grey90",border=F)
for(i  in 1:(nstrat))
{
  if(i==8) next
  h<-i
  if(i==9) h<-8
  points(0.05,h,pch=16,cex=1.2,col=cols[i])
  lines(c(0.02,.08),c(h,h),col=cols[i],lwd=1)
  text(.08,h,stratlabs[i],cex=1.2,pos=4)
  
}
text(0.05,8.5,"(b)")
dev.off()
