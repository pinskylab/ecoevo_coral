#==========================================================================
#
# Create heatmaps for coral management strategies
#
# - genetic variance (V) on the x-axis
# - dispersal rate (D) on the y-axis
# - one panel per management strategy
# - Values at each combination of D and V are coral cover relative to
#     the coral cover under the best performing strategy for each scenario
#==========================================================================


# Read in data
#load("C:/Users/Tim Walsworth/Dropbox/Corals2/SA_RandomTemp_testinng_CoralRunsOut031017.RData") # located in Corals2 folder
iterout<-read.csv("Base_HeatmapsReserveSizeSimulations_AverageDensities_newStrats_2017-05-16.csv",header=T)
iterout2<-read.csv("Base_HeatmapsReserveSizeSimulations_RelativeDensities_newStrats_2017-05-16.csv",header=T)
head(scenarios)

stratlabs<-c("Hot","Cold","Hot and Cold","High Cover","Low Cover","Gradient","Portfolio","Portfolio 2","Random")

strat1<-iterout2[iterout2[,2]==1,]
res20<-iterout2[iterout2[,3]==0.2,]
av20<-iterout[iterout[,3]==.2,]
xtab1<-xtabs(AvDensity~V+D+Strategy,data=res20)/50
xtab12<-xtabs(AvDensity~V+D+Strategy,data=av20)/50
xtab2<-xtab1
xtab3<-xtab2
for(i in 1:9)
{
  for(j in 1:3)
  {
    for(k in 1:3)
    {
      xtab2[j,k,i]<-(xtab1[j,k,i]-min(xtab1[j,k,]))/(max(xtab1[j,k,]-min(xtab1[j,k,])))
      #xtab3[j,k,i]<-1-(max(xtab12[j,k,])-xtab12[j,k,i])
      xtab3[j,k,i]<-(xtab12[j,k,i]/max(xtab12))
    }
  }
  
}
#xtab3<-1-((1-xtab3)/max(1-xtab3))

cols<-colorRampPalette(c("darkblue","red"))(100)
breakpts<-seq(0,1,by=.01)
#cols<-colorRampPalette(c("white","darkred"))(100)

colbreaks<-seq(.01,1,by=.01)
boxtype<-c("n","u","n","]","n","c","7","7","n")

pdf("SafeOperatingSpace_Rainbow.pdf",width=10,height=8,pointsize=12)
layout(matrix(c(1,1,10,10,2,2,11,11,3,3,
                1,1,10,10,2,2,11,11,3,3,
                12,12,12,12,12,12,12,12,12,12,
                4,4,10,10,9,9,11,11,5,5,
                4,4,10,10,9,9,11,11,5,5,
                13,13,13,13,13,13,13,13,13,13,
                6,6,10,10,7,7,11,11,8,8,
                6,6,10,10,7,7,11,11,8,8),nrow=8,ncol=10,byrow=T),heights=c(1,1,.05,1,1,.05,1,1.2),
       widths=c(1,1,.05,.05,1,1,.05,.05,1,1))
par(oma=c(2,4,0,4),xpd=TRUE,fg="black",cex.axis=1.5,cex.lab=3)
for(i in c(6,7,9,4,5,1,2,3))
{
  if(i==8) next
  par(mar=c(1,1,1.8,1))
  if(i<4) par(mar=c(3,1,1.8,1))
  # strat1<-iterout2[iterout2[,2]==i,]
  # xtab1<-xtabs(AvDensity~V+D,data=strat1)
  # print(xtab1)
  # # a<-image(1-xtab1,col=brewer.pal(11,"Spectral"),axes=F)
  # print(a)
  # if(i == 7)
  # {
  #   axis(side=1,at=seq(0,1,length.out=ncol(xtab1)),labels=colnames(xtab1))
  #   axis(side=2,at=seq(1,0,length.out=nrow(xtab1)),labels=rownames(xtab1))
  #   mtext(side=1,line=3,"D")
  #   mtext(side=2,line=3,"V")
  # }
  #filled.contour3(xtab2[,,i],col=cols,zlim=c(0,1))
  #image(xtab2[,,i],col=cols,axes=F,main=strategies[i],breaks=breakpts)
  plot(NA,NA,xlim=c(0,3),ylim=c(0,3),
       yaxs="i",xaxs="i",bty="o",axes=F,xpd=T,main=stratlabs[i],ylab="",xlab="",cex.main=1.5)
  for(q in 1:3)
  {
    for(v in 1:3)
    {
      # rad<-sqrt(xtab3[q,v,i]/pi)
      
      inten<-xtab3[q,v,i]
      
      xvals<-c(q-1,q,q,q-1)
      yvals<-c(v-1,v-1,v,v)
      
      # theta<-seq(0,2*pi,length.out=1000)
      # x.center<-sqrt(1/pi)+(q-1)*2*sqrt(1/pi)
      # y.center<-sqrt(1/pi)+(v-1)*2*sqrt(1/pi)
      colindex<-which.min(abs(colbreaks-xtab2[q,v,i]))
      print(colindex)
      if(colbreaks[colindex]-xtab2[q,v,i]<0)colindex<-colindex+1
      polygon(x=xvals,y=yvals,col=adjustcolor(cols[colindex],alpha.f=inten),
              border="grey",lwd=.8)
      
      
    }
  }
  #image(xtab2[,,i],col=cols,axes=F,main=strategies[i])
  if(i<4)  axis(side=1,at=c(0,1,2)+.5,labels=c(0,.1,.4),cex.axis=1.5)
  if(i %in% c(1,4,6)) axis(side=2,at=c(0,1,2)+.5,labels=c(0,.001,.01),cex.axis=1.5)
  labs<-as.vector(signif(xtab2[,,i],2))
  labs[labs<0.01]<-0
 # text(y=c(0,0,0,.5,.5,.5,1,1,1),x=rep(c(0,.5,1),3),labels=labs)
  #text(y=c(1,1,1,.5,.5,.5,0,0,0),x=rep(c(0,.5,1),3),labels=seq(1,9))
  #box("figure",col="grey",bty=boxtype[i])
  # if(i==3) text(sqrt(1/pi)+(4-1)*2*sqrt(1/pi),sqrt(1/pi)+(2-1)*2*sqrt(1/pi),"Temperature",srt=270,cex=1.5)
  # if(i==5) text(sqrt(1/pi)+(4-1)*2*sqrt(1/pi),sqrt(1/pi)+(2-1)*2*sqrt(1/pi),"Cover",srt=270,cex=1.5)
  # if(i==9) text(sqrt(1/pi)+(4-1)*2*sqrt(1/pi),sqrt(1/pi)+(2-1)*2*sqrt(1/pi),"Diversity",srt=270,cex=1.5)
  }
mtext(side=2,outer=T,line=2,"Dispersal",cex=1.5)
mtext(side=1,outer=T,line=1,"Genetic Variance",cex=1.5)
xval<-seq(0,100,by=1)
par(mar=c(4,4,1,1))
plot(seq(0,100),seq(0,100),type="n",ylab="",xlab="",yaxt="n",xaxt="n",
     xaxs="i",yaxs="i",bty="o",ylim=c(0,100),xlim=c(0,100))
for(i in 1:100)
{
  for(z in 1:100)
  {
    polygon(c(i-1,i,i,i-1),c(z-1,z-1,z,z),col=adjustcolor(cols[i],alpha.f=z/100),lty="blank",border=F)
  }
  
}
axis(side=1,at=c(0,100),labels=c("Worst","Best"))
axis(side=2,at=c(0,100),labels=c(0,1))
mtext(side=2,line=3,"Relative Coral Cover",cex=1.2)
mtext(side=1,line=3,"Relative Management Performance",cex=1.2)

par(mar=c(0,0,0,0),fg="grey")
plot(1,1,type="n",axes=F)
abline(v=1)
plot(1,1,type="n",axes=F)
abline(v=1)
plot(1,1,type="n",axes=F)
abline(h=1)
plot(1,1,type="n",axes=F)
abline(h=1)
mtext(side=4,outer=T,line=0.2,cex=1.5,at=rev(c(.84,.5,.18)),c("Temperature","Cover","Diverse Portfolio"),col="black")
mtext(side=4,outer=T,line=1.8,cex=1.5,at=c(.84,.5,.18),"Strategies",col="black")
dev.off()

breakpts2<-seq(min(xtab3),0,length.out=101)
pdf("SafeOperatingSpace_diffFromBest.pdf",width=8,height=6)
layout(matrix(c(10,10,10,
                1,2,3,
                4,5,6,
                7,8,9),nrow=4,ncol=3,byrow=T),heights=c(0.5,1,1,1))
par(oma=c(2,4,0,4))
for(i in 1:9)
{
  par(mar=c(1,1,1.5,1))
  if(i>6) par(mar=c(3,1,1.5,1))
  # strat1<-iterout2[iterout2[,2]==i,]
  # xtab1<-xtabs(AvDensity~V+D,data=strat1)
  # print(xtab1)
  # # a<-image(1-xtab1,col=brewer.pal(11,"Spectral"),axes=F)
  # print(a)
  # if(i == 7)
  # {
  #   axis(side=1,at=seq(0,1,length.out=ncol(xtab1)),labels=colnames(xtab1))
  #   axis(side=2,at=seq(1,0,length.out=nrow(xtab1)),labels=rownames(xtab1))
  #   mtext(side=1,line=3,"D")
  #   mtext(side=2,line=3,"V")
  # }
  #filled.contour3(xtab2[,,i],col=cols,zlim=c(0,1))
  image(xtab3[,,i],col=cols,axes=F,main=strategies[i],breaks=breakpts2)
  #image(xtab2[,,i],col=cols,axes=F,main=strategies[i])
  if(i>6)  axis(side=1,at=c(0,.5,1),labels=c(0,.1,.4),cex.axis=1.2)
  if(i %in% c(1,4,7)) axis(side=2,at=c(0,.5,1),labels=c(0,.001,.01),cex.axis=1.2)
  labs<-as.vector(signif(xtab3[,,i],3))
  labs[abs(labs)<0.01]<-0
  text(y=c(0,0,0,.5,.5,.5,1,1,1),x=rep(c(0,.5,1),3),labels=labs)
}
mtext(side=2,outer=T,line=2,"Dispersal",cex=1.5)
mtext(side=1,outer=T,line=1,"Genetic Variance",cex=1.5)
xval<-seq(0,100,by=1)
par(mar=c(1,1,1.5,1))
plot(xval,rep(0,101),type="n",yaxt="n",xaxt="n",bty="n",ylim=c(-1,1))
for(i in 1:100)
{
  polygon(c(xval[i],xval[i+1],xval[i+1],xval[i]),c(0,0,1,1),col=cols[i],border=F)
}
text(5,-.5,"Minimum",cex=1.2)
text(95,-.5,"Maximum",cex=1.2)
text(50,-.5,"Relative Performance",cex=1.5)


dev.off()


cols[1]<-"black"
pdf("SafeOperatingSpace_RelativetoBestandWorst_blackout.pdf",width=8,height=6)
layout(matrix(c(10,10,10,
                1,2,3,
                4,5,6,
                7,8,9),nrow=4,ncol=3,byrow=T),heights=c(0.5,1,1,1))
par(oma=c(2,4,0,4))
for(i in 1:9)
{
  par(mar=c(1,1,1.5,1))
  if(i>6) par(mar=c(3,1,1.5,1))
  # strat1<-iterout2[iterout2[,2]==i,]
  # xtab1<-xtabs(AvDensity~V+D,data=strat1)
  # print(xtab1)
  # # a<-image(1-xtab1,col=brewer.pal(11,"Spectral"),axes=F)
  # print(a)
  # if(i == 7)
  # {
  #   axis(side=1,at=seq(0,1,length.out=ncol(xtab1)),labels=colnames(xtab1))
  #   axis(side=2,at=seq(1,0,length.out=nrow(xtab1)),labels=rownames(xtab1))
  #   mtext(side=1,line=3,"D")
  #   mtext(side=2,line=3,"V")
  # }
  #filled.contour3(xtab2[,,i],col=cols,zlim=c(0,1))
  image(xtab2[,,i],col=cols,axes=F,main=strategies[i],breaks=breakpts)
  #image(xtab2[,,i],col=cols,axes=F,main=strategies[i])
  if(i>6)  axis(side=1,at=c(0,.5,1),labels=c(0,.1,.4),cex.axis=1.2)
  if(i %in% c(1,4,7)) axis(side=2,at=c(0,.5,1),labels=c(0,.001,.01),cex.axis=1.2)
  labs<-as.vector(signif(xtab2[,,i],2))
  labs[labs<0.01]<-0
  text(y=c(0,0,0,.5,.5,.5,1,1,1),x=rep(c(0,.5,1),3),labels=labs)
  #text(y=c(1,1,1,.5,.5,.5,0,0,0),x=rep(c(0,.5,1),3),labels=seq(1,9))
}
mtext(side=2,outer=T,line=2,"Dispersal",cex=1.5)
mtext(side=1,outer=T,line=1,"Genetic Variance",cex=1.5)
xval<-seq(0,100,by=1)
par(mar=c(1,1,1.5,1))
plot(xval,rep(0,101),type="n",yaxt="n",xaxt="n",bty="n",ylim=c(-1,1))
for(i in 1:100)
{
  polygon(c(xval[i],xval[i+1],xval[i+1],xval[i]),c(0,0,1,1),col=cols[i],border=F)
}
text(5,-.5,"Minimum",cex=1.2)
text(95,-.5,"Maximum",cex=1.2)
text(50,-.5,"Relative Performance",cex=1.5)


dev.off()


my.filled.contour<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                          length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                             ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                             levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                             col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                             key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                             axes = TRUE, frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  #===========================================================
  
  #   layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  #   plot.new()
  #   plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
  #               yaxs = "i")
  #   rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  #   if (missing(key.axes)) {
  #     if (axes) 
  #       axis(4)
  #   }
  #   else key.axes
  #   box()
  #   if (!missing(key.title)) 
  #     key.title
  #   #===========================================================
  
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


diverge.color <- function(data,pal_choice="RdBu",centeredOn=0){
  nHalf=50
  Min <- min(data,na.rm=TRUE)
  Max <- max(data,na.rm=TRUE)
  Thresh <- centeredOn
  pal<-brewer.pal(n=11,pal_choice)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:10){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- pretty(Min, Thresh, length.out=nHalf+1)
  rb2 <- pretty(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  cuts <- classIntervals(data, style="fixed",fixedBreaks=rampbreaks)
  return(list(cuts,rc1))
}

filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }


