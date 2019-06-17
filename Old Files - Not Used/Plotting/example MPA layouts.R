# Making an example MPA layout figure

load("C:/Users/Tim Walsworth/Dropbox/AWSCoralRuns/Base_ReserveSize_NewStrategiesExtinction_CoralSimulationsMPASetOutData012317.RData")

stratlabs<-c("Hot","Cold","Hot and Cold","High Cover","Low Cover","Gradient","Portfolio","Portfolio 2","Random")

oneiter<-mpaset[,,1,]

onescen<-oneiter[,,26]
onescen<-onescen[,-8]

layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(1,8))
par(mar=c(0,0,0,0))
frame()
legend("center",pch=15,legend=c("No MPA","MPA"),col=c("grey","dodgerblue"),
       horiz = T,cex=1.5,bty="n",pt.cex=3)
par(mar=c(5,8,0,2))
image(onescen[,seq(8,1)],col=colorRampPalette(c("grey","dodgerblue"))(10),axes=F)
axis(side=1,at=c(0.5,19.5,39.5,59.5)/60,labels=c(1,20,40,60))
axis(side=2,at=seq(0,1,length.out=8),labels=rev(stratlabs[-8]),las=2,cex.axis=1.2)
mtext(side=1,line=3,"Location",cex=1.5)
