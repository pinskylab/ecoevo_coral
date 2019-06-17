#==========================================================================
#
# Create figure of coral cover inside and outside of reserves under different
#  scenarios of evolution and thermal landscape arrangement. This is figure
#  2 in Walsworth et al. 2019 Nature Climate Change
#
#==========================================================================
load("Walsworthetal_NCC_Figure2Data.Rdata")

cols1<-c("lightblue","dodgerblue","dodgerblue4")
cols2<-colorRampPalette(c("burlywood","darkorange"))(2)
cols3<-colorRampPalette(c("green2","darkgreen"))(3)[c(3,1,2)]
cols<-c(cols1,cols2,cols3)

stratlabs<-c("Hot Reefs","Cold Reefs","Hot and Cold\nReefs","High Cover","Low Cover","Evenly-spaced","Portfolio","Random")

pdf("Walsworthetal_Fig2_Proof.pdf",height=10,width=9)
layout(matrix(c(1,2,
                3,4,
                5,5),nrow=3,ncol=2,byrow=T),heights=c(1,1,.3))
par(xpd=F,mar=c(5,5,1,1),cex.lab=1.3,cex.text=1.8,cex.axis=1.2,cex.legend=1.5)

#===========================
# No genetic variance
#===========================
plot(NA,NA,xlim=c(0,.8),ylim=c(0,.8),xaxs="i",yaxs="i",xlab="Proportion Coral Cover Inside Reserves",
     ylab="Proportion Coral Cover Outside Reserves")
polygon(c(0,1,1,0),c(0,0,1,0),col="grey80",border=F)

for(i in 1:8)
{
  lines(incovall1[,i],outcovall1[,i],col=cols[i],lwd=2)
}

text(.3,.35,"0",cex=1.5)
text(.57,.07,"100-500",cex=1.5)
segments(0,0,.51,.07)
legend("topleft",bty="n",legend="(a) No Genetic Variance",cex=1.5)

#========================
# Low genetic variance
#========================
plot(NA,NA,xlim=c(0,.8),ylim=c(0,.8),xaxs="i",yaxs="i",xlab="Proportion Coral Cover Inside Reserves",
     ylab=" Proportion Coral Cover Outside Reserves")
polygon(c(0,1,1,0),c(0,0,1,0),col="grey80",border=F)


for(i in 1:8)
{
  lines(incovall2[,i],outcovall2[,i],col=cols[i],lwd=2)
}

text(rowMeans(incovall2[c(1,100,200,300,400,500),])+c(-.03,-.05,-.05,0,0,.03),
     rowMeans(outcovall2[c(1,100,200,300,400,500),])+c(0,0,.02,.05,.05,.0),
       labels=c("0","100","200","300","400","500"),cex=1.5)

legend("topleft",bty="n",legend="(b) Low Genetic Variance",cex=1.5)

#==============================
# High Genetic Variance
#==============================
plot(NA,NA,xlim=c(0,.8),ylim=c(0,.8),xaxs="i",yaxs="i",xlab="Proportion Coral Cover Inside Reserves",
     ylab="Proportion Coral Cover Outside Reserves")
polygon(c(0,1,1,0),c(0,0,1,0),col="grey80",border=F)

for(i in 1:8)
{
  lines(incovall3[,i],outcovall3[,i],col=cols[i],lwd=2)
}

text(rowMeans(incovall3[c(1,100,200,300,400,500),])+c(0,-.05,-.02,0,0.08,.04),
     rowMeans(outcovall3[c(1,100,200,300,400,500),])+c(0,0,-.02,-.05,0,.0),
     labels=c("0","100","200","300","400","500"),cex=1.5)

legend("topleft",bty="n",legend="(c) High Genetic Variance",cex=1.5)

#===============================================================
# Randomized thermal landscape and low genetic variance
#===============================================================
plot(NA,NA,xlim=c(0,.8),ylim=c(0,.8),xaxs="i",yaxs="i",xlab="Proportion Coral Cover Inside Reserves",
     ylab="Proportion Coral Cover Outside Reserves")
polygon(c(0,1,1,0),c(0,0,1,0),col="grey80",border=F)


for(i in 1:8)
{
  lines(incovall4[,i],outcovall4[,i],col=cols[i],lwd=2)
}

text(rowMeans(incovall4[c(1,100,200,300,400,500),])+c(0,0.02,-.04,0,0.02,0),
     rowMeans(outcovall4[c(1,100,200,300,400,500),])+c(0.02,0,-.04,-.02,-.03,.02),
     labels=c("0","100","200","300","400","500"),cex=1.5)


legend("topleft",bty="n",legend="(d) Randomized Thermal Landscape\nLow Genetic Variance",cex=1.5)

par(mar=c(1,5,1,1))
frame()
legend("center",col=cols,legend=stratlabs,cex=1.5,lwd=2,lty=1,bty="n",ncol=4)

dev.off()
