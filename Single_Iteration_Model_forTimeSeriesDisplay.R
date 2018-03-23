source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1)
scenarios<-read.csv("Norberg_StartingParameters_iterationlist_reservesize.csv",header=T)
#start.parms<-start.parm.mat
start.parms<-matrix(scenarios[8,-c(1,17,18)],nrow=5,ncol=3,byrow=F)
colnames(start.parms)<-colnames(start.parm.mat)
rownames(start.parms)<-rownames(start.parm.mat)

library(mvtnorm)


# Making A Change
set.seed(6)
nsp<-3                        # how many species in model?
size<-60                      # how many reefs in model?
maxtime1<-1500 
burnin<-maxtime1
runtime<-500# how many time steps?
times<-seq(0,maxtime1,by=1)    #Vector from 1 to the number of time steps
mid<-27                       # mean temperature across all reefs at start of simulation.
range<-3                      # range of temperatures across reefs at start of simulation
MPAamount<-0.2
cgrow<-1
agrow<-1
tempchange<-3
monitor.yrs<-c(50,100,500)

temp.stoch<-1.2
mdim<-size
lindec<-exp(seq(0,-5,length=mdim))#rev(seq(-0.9,1,length=mdim))
ma<-matrix(0,nrow=mdim,ncol=mdim)
ma[,1]<-lindec
for(i in 2:mdim){
  ma[-c(1:(i-1)),i]<-lindec[1:(mdim-(i-1))]
}
ma<-ma+t(ma)-diag(nrow(ma))
is.positive.definite(ma)


sds<-rep(.2*temp.stoch,size)
b<-sds%*%t(sds)
spatialtemp<-b*ma

anoms.burn<-matrix(rnorm(burnin,0,temp.stoch),nrow=size,ncol=burnin,byrow=T)+matrix(rmvnorm(burnin,rep(0,size),spatialtemp),nrow=size,byrow=T)

anoms.runs<-matrix(rnorm(runtime,0,temp.stoch),nrow=size,ncol=runtime,byrow=T)+matrix(rmvnorm(runtime,rep(0,size),spatialtemp),nrow=size,byrow=T)

algaemort<-matrix(runif((runtime+burnin)*size,0.05,.3),nrow=size,ncol=runtime+burnin,byrow=T)

#==============================================================================
# Generate starting conditions and merge into single vector or use in ode.1D()
#   - Temperatures (specify scenario)
#   - Densities by species
#   - Trait (optimum temperatures, specify scenario)
#
# Create vector of names for each "State"
#   - species names, then trait names, then temperature
#==============================================================================
temps<-generate.temps(size=size,mid=mid,range=range,temp.scenario="linear") 
sppstate<-generate.state(size=size,nsp=nsp,dens=0.25,random=T)
traitstate<-generate.traits(nsp=nsp,size=size,mid=mid,temps=temps,range=range,trait.scenario="perfect.adapt")
allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)

allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")

#==============================================================================
# What types of species (will they be protected by an MPA?)
#   - Set species that are protected by MPA to sptype=1
#   - Set species that are NOT protected by MPA to sptype=2
#   - For species where sptype=2, set mpa to a value >1
#   - For species where sptype=1, set mpa<=1
#   - mpa will be multiplied by mortality rate
#==============================================================================
species<-c("C1","C2","MA")
sptype<-c(1,1,2)
mpa<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy="none")
# mpa<-rep(0,size)
# mpa[c((round(size/3)):(2*round(size/3)))]<-1

#==============================================================================
# Define model parameters (# required)
#   - nsp: Number of Species (1)
#   - V: Genetic variance (nsp)
#   - D: Dispersal (nsp)
#   - rmax: Maximum growth rate (nsp)
#   - alphas: species interaction (competition) matrix (nsp x nsp matrix)
#   - m.normal: mortality rate in normal year (nsp)
#   - m.catastrophe: mortality rate in catastrophic year (nsp)
#   - pcatastrophe: annual probability of catastrophic mortality (1)
#   - w: temperature tolerance or plasticity (nsp)
#   - annual.temp.change: rate of temperature change (1)
#   - maxtemp: maximum temperature (sigmoid temp change scenario) (1)
#   - Nmin: minimum density (for computational happiness) (1)
#   - deltax: step size in space for derivitives (1)
#==============================================================================

parms<-list(
  nsp=nsp,
  #mpa=mpa,
  monitor.yrs<-monitor.yrs,
  species=species,
  sptype=sptype,
  V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
  D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
  rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
  #alphas=matrix(1,nrow=nsp,ncol=nsp),
  #alphas=diag(1,nrow=nsp,ncol=nsp),
  alphas=matrix(c(1,1.1,1,
                  1,1,1,
                  .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
  m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
  
  w=c(t(as.numeric(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
  pcatastrophe=0.02,
  annual.temp.change=.011,
  maxtemp=30,
  Nmin=10^-6,
  deltax=1,
  spatialtemp=spatialtemp,
  timemod=0
)

#===================================================================================
# Function to calculate partial derivitives of population growth, trait change, and 
#     temperature change at a given time step. Appends the new state values to 
#     a matrix storing a time serie sof state values. Additionally creates a time 
#     series of state change values. Returns a list with the two time series
#     matrices.
#
# Parameters:
#     t: time step
#     y: Vector of state values (input "allstate" vector here)
#     parms: List of parameters to be used in the model
#     size: Size of reef
#     nsp: Number of species
#     temp.change: Scenario of temperature change to be modeled. "const" has
#                   constant temperatures through time. "linear" has linear
#                   temperature change. "sigmoid" has logistic temperature change up
#                   to a maxtemp value.
#====================================================================================

coral_trait_stoch<-function(t,y, parms,size,nsp,temp.change=c("const","linear","sigmoid"),
                            stoch.temp=TRUE,anoms,algaemort){
  
  with(as.list(parms), {  # Load parameter values
    print("A")
    time.series<-matrix(NA,nrow=(2*nsp*size+size),ncol=t) # Create storage matrix for time series of state values at each location
    dspout<-time.series # Create storage matrix for time series of change in state values at each location
    time.series[,1]<-y  # set initial state
    extinctions<-matrix(0,nrow=3,ncol=(t-1))
    for(k in 2:t)
    {
      print(k) # Print counter
      
      spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
      traits<-matrix(time.series[(size*nsp+1):(length(y)-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
      temps<-time.series[(length(y)-(size-1)):length(y),k-1]+anoms[,k-1] # Extract temperature state values
      #print(temps)  
#       pcat<-rbinom(1,1,pcatastrophe)           # Binomial draw for annual catastrophe occurance
#       m<-(1-pcat)*m.normal+pcat*m.catastrophe  # set annual mortality rate dependent on presence of catastrophe
      print("B")  
      m<-m
      dspp<-spps   # create storage for change in density state values
      dtraits<-traits # create storage for change in trait state values
      #print(w)
      for(i in 1:nsp)
      { 
        print("C")
        #print(w[i]^2)
        if(i==3) amort<-algaemort[,k+timemod]
        dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                       mort=amort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
        dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                          mort=amort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                          mortality.model="tempvary",spp=species[i],growth="normal")
      }
      
      # Calculate change in temperatures across reef
      dtemps<-switch(temp.change,
                     const=rep(0,size),
                     linear=rep(annual.temp.change,size),
                     sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))
      
      dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector
      

      dspout[,k]<-dsp  # Store change in state value
      time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
      # extinctions[,k-1]<-extinctrisk(Nall=matrix(time.series[1:(3*size),k],nrow=3,ncol=size,byrow=T),
      #                                threshold=.01)
    }
    return(list("ts"=time.series,"dsp"=dspout,"ext"=extinctions))
  })


  
}
start.time<-Sys.time()
#out<-apply(allstate,MARGIN=2,FUN=coral_trait_stoch,t=maxtime1,parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
 #           temp.change="const",stoch.temp=T,burn=T,anoms=anoms.burn,amort=amort)
out<-coral_trait_stoch(t=maxtime1,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                       temp.change="const",stoch.temp=T,anoms=anoms.burn,algaemort=algaemort)
Sys.time()-start.time

#   }
# }

temps<-out$ts[(2*nsp*size+1):(2*nsp*size+size),maxtime1] 
sppstate<-matrix(out$ts[1:(nsp*size),maxtime1],nrow=nsp,ncol=size,byrow=T) 
traitstate<-out$ts[(nsp*size+1):(2*nsp*size),maxtime1] 
allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)

allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")
cols<-colorRampPalette(c("red","yellow","blue"))
cols2<-c("dodgerblue","darkorange")
strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","random") # Name different management strategies
#strategies<-c("portfolioHC","portfolio","random","none")
nstrat<-length(strategies)   
mpa<-matrix(NA,nrow=size,ncol=length(strategies))
out2<-array(NA,dim=c(420,500,length(strategies)))
extinctionout<-array(NA,dim=c(3,499,length(strategies)))
for(i in 1:length(strategies))
{
  if(strategies[i]!="switch"){
    mpa[,i]<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy=strategies[i],priordata=out$ts)
    maxtime2<-500
    
  }
  print(mpa)
  # if(strategies[i]=="switch")
  # {
  #   mpa[,i]<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy=strategies[2],priordata=out$ts)
  #   maxtime2<-200
  # }
  parms<-list(
    nsp=nsp,
    mpa=mpa[,i],
    monitor.yrs<-monitor.yrs,
    species=species,
    sptype=sptype,
    V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
    D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
    rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
    #alphas=matrix(1,nrow=nsp,ncol=nsp),
    #alphas=diag(1,nrow=nsp,ncol=nsp),
    alphas=matrix(c(1,1.1,1,
                    1,1,1,
                    .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
    m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
    
    w=c(t(as.numeric(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
    pcatastrophe=0.02,
    annual.temp.change=.011,
    maxtemp=30,
    Nmin=10^-6,
    deltax=1,
    spatialtemp=spatialtemp,
    timemod=burnin
  )
  
  modout<-coral_trait_stoch(t=maxtime2,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                              temp.change="sigmoid",stoch.temp=T,anoms=anoms.runs,algaemort=algaemort)
  out2[,1:maxtime2,i]<-modout$ts
  extinctionout[,1:(maxtime2-1),i]<-modout$ext
  
  # if(strategies[i]=="switch")
  # {
  #   temps<-out2[(2*nsp*size+1):(2*nsp*size+size),maxtime2,i] 
  #   sppstate<-matrix(out2[1:(nsp*size),maxtime2,i],nrow=nsp,ncol=size,byrow=T) 
  #   traitstate<-out2[(nsp*size+1):(2*nsp*size),maxtime2,i] 
  #   allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)
  #   mpa[,i]<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy="portfolioGreedy",priordata=out$ts)
  #   maxtime3<-300
  #   parms<-list(
  #     nsp=nsp,
  #     mpa=mpa[,i],
  #     V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
  #     D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
  #     rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
  #     #  alphas=matrix(1,nrow=nsp,ncol=nsp),
  #     #alphas=diag(1,nrow=nsp,ncol=nsp),
  #     alphas=matrix(c(1,1.1,1,
  #                     1,1,1,
  #                     .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
  #     m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
  #     
  #     w=c(t(as.numeric(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
  #     pcatastrophe=0.02,
  #     temp.stoch=1.2, #0.9 works great for constant temp
  #     annual.temp.change=.011,
  #     maxtemp=mid+tempchange,
  #     Nmin=10^-6,
  #     deltax=1,
  #     spatialtemp=spatialtemp
  #     
  #   )
  #   modout<-coral_trait_stoch(t=maxtime3,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
  #                             temp.change="sigmoid",stoch.temp=T,anoms=anoms.runs[,-c(1:maxtime2)])
  #   out2[,(maxtime2+1):(maxtime2+maxtime3),i]<-modout$ts
  #   dim(modout$ext)
  #   extinctionout[,(maxtime2):(maxtime2+maxtime3-2),i]<-modout$ext
  # }
  
#   j<-1 # start counter
#   while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
#   {
# 
#     lines(seq(1,runtime)+1500,(colSums(out$ts[((j-1)*size+1):(j*size),])/size),col=cols2[j],lwd=1)
#     j<-j+1
#   }
# lines(seq(1,runtime)+1500,(colSums(out$ts[((1-1)*size+1):(max(hard.corals)*size),])/size),col=cols(7)[i],lwd=2) # Plot the total hard coral density across the reef through time

}
# 
# 
# library(gtools)
# 
# combinations(60,20,v=1:60,set=T,repeats.allowed = F)
# 
# 
# 
# 
# 
# 
# 
# 
# ti<-Sys.time()
# aa<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy="portfolioHC",priordata=out$ts)
# Sys.time()-ti
# 
# a<-array(sppstate,dim=c(3,size,50))
# 
# ti<-Sys.time()
# apply(a,MARGIN=3,FUN=setMPA,temps=temps,sptype=sptype,size=size,amount=MPAamount,strategy="portfolio",priordata=out$ts)
# Sys.time()-ti
# 
# setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy=strategies[i],priordata=out$ts)
# 
# dzbar<-array(NA,dim=c(60,500,length(strategies)))
# dzbar2<-dzbar
# for(z in 1:length(strategies))
# {
#   for(i in 1:60)
#   {
#     eco<-dZbardt(Nall=out2[c((0+i),(60+i)),,z],Zall=out2[c((180+i),(240+i)),,z])$Ecology
#     evo<-dZbardt(Nall=out2[c((0+i),(60+i)),,z],Zall=out2[c((180+i),(240+i)),,z])$Evolution
#     dzbar[i,,z]<-(evo/(eco+evo))
#     dzbar2[i,,z]<-eco+evo
#   }
#   
# }
# 
# colramp<-colorRampPalette(c("dodgerblue","white","darkorange"))(50)
# 
# layout(matrix(c(1,1,
#                 2,3,
#                 4,5,
#                 6,7,
#                 8,9),nrow=5,ncol=2,byrow=T),heights=c(1,4,4,4,4))
# par(mar=c(1,3,1,1))
# plot(rep(1,50),col=colramp,pch=16,xaxt="n",yaxt="n",bty="n",cex=3)
# axis(side=1,at=c(1,50),labels=c("Eco","Evo"))
# par(mar=c(4,5,1,2)+.1)
# for(i in 1:length(strategies))
# {
#   image(t(dzbar[,,i]),col=colramp,xaxt="n",yaxt="n",ylab="Location",xlab="Time")
#   axis(side=1,at=seq(0,1,length.out=6),labels=seq(0,500,length.out=6))
#   axis(side=2,at=seq(0,1,length.out=6),labels=seq(0,60,length.out=6))
#   text(.1,.9,labels=strategies[i],cex=2)
# }
# 
# image(t(dzbar[,,2]),col=adjustcolor(colramp,
#                                     alpha.f=t(dzbar2[,,2]/max(abs(dzbar2[,,2])))),xlim=c(0,500))
# 
# a<-dZbardt(Nall=out2[c(1,61),,5],Zall=out2[c(181,241),,5])
# 
# plot(a$Ecology,type="l",lwd=2,col="darkgreen")
# lines(a$Evolution,lwd=2,col="darkorange")
# 
# 
# 
# sp1cover<-out2[1:60,,]
# sp2cover<-out2[61:120,,]
# sp1trait<-out2[181:240,,]
# sp2trait<-out2[241:300,,]
# wtd.var(x=sp1trait[,500,7],weights=sp1cover[,500,7])
# wtd.mean(x=sp1trait[,500,7],weights=sp1cover[,500,7])
# 
# wtstats<-array(NA,dim=c(length(strategies),3,500))
# 
# for(i in 1:length(strategies))
# {
#   for(j in 1:500)
#   {
#     wtstats[i,1,j]<-wtd.mean(sp1trait[,j,i],weights=sp1cover[,j,i],normwt=T)
#     wtstats[i,2,j]<-wtd.var(sp1trait[,j,i],weights=sp1cover[,j,i],normwt=T)
#     wtstats[i,3,j]<-mean(sp1cover[,j,i])
#   }
# }
# 
# library(RColorBrewer)
# cols<-c(brewer.pal(7,"Dark2")[c(2,1,3,4,5,6,7)],"black")
# layout(matrix(c(1,1,2,2,
#                 1,1,2,2,
#                 1,1,2,2,
#                 1,1,2,2,
#                 3,3,4,4,
#                 3,3,4,4,
#                 3,3,4,4,
#                 3,3,4,4),nrow=8,ncol=4,byrow=T))
# plot(1,1,type="n",xlim=c(26,36),ylim=c(0,20),bty="l",ylab="Weighted Variance Trait",
#      xlab="Weighted Mean Trait")
# 
# for(i in c(8,1,2,3,4,5,6,7))
# {
#   tempcols<-colorRampPalette(c("grey",cols[i]))(500)
#   afs<-wtstats[i,3,]/max(wtstats[,3,])
#  # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16,cex=seq(1,500,by=1)/125)
#   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16)
#   for(j in 2:500)
#   {
#     points(wtstats[i,1,j],wtstats[i,2,j],
#            col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),
#            pch=16,cex=j/125)
#     arrows(wtstats[i,1,j-1],wtstats[i,2,j-1],wtstats[i,1,j],wtstats[i,2,j],lwd=4,
#            col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),length=0.05)
#   }
# }
# 
# 
# for(i in c(7,1,2,3,4,5,6))
# {
#   print(i)
#   if(i==7) plot(density(diff(wtstats[i,2,])),col=adjustcolor(cols[i],alpha.f=.8),lwd=2,
#                 bty="n",main="")
#   else lines(density(diff(wtstats[i,2,])),col=adjustcolor(cols[i],alpha.f=.8),lwd=2)
#     
# }
# 
# #legend("topleft",horiz=F,legend=strategies,pch=16,col=cols,bty="n",cex=2)
# 
# # par(mar=c(4,4,0,0))
# # for(i in 1:7)
# # {
# #   plot(1,1,type="n",xlim=c(1,500),ylim=c(-.4,.4),bty="l",ylab=expression(paste(Delta, "Weighted Variance Trait")),xlab="Time")
# #   aa<-diff(wtstats[i,2,])
# #   tempcols<-colorRampPalette(c("grey",cols[i]))(500)
# #   #afs<-wtstats[i,3,]/max(wtstats[,3,])
# #   lines(seq(1,499),aa,col=tempcols[500],lwd=1)
# #   polypoints<-rep(0,499*2)
# #   
# #   for(q in 1:(499))
# #   {
# #     polypoints[q]<-max(aa[q],0)
# #   }
# #   
# #   for(q in (2*499):500)
# #   {
# #     print(q)
# #     polypoints[q]<-min(aa[q-499],0)
# #   }
# #   bb<-cbind(c(seq(1,499,by=1),seq(499,1,by=-1)),polypoints)
# #   polygon(bb, col=adjustcolor(tempcols[500],alpha.f=.8),border=F)
# #   
# #   abline(h=0,lwd=.5)
# #   abline(h=mean(diff(wtstats[i,2,])),lty=2,lwd=2)
# #   #print(mean(diff(wtstats[i,2,])))
# #   abline(h=quantile(diff(wtstats[i,2,]),c(.025,.975))[c(1,2)])
# #   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=afs),pch=16,cex=seq(1,500,by=1)/125)
# #   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16)
# #   # for(j in 2:500)
# #   # {
# #   #   
# #   #   # points(wtstats[i,1,j],wtstats[i,2,j],col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),pch=16,cex=wtstats[i,3,j]/max(wtstats[,3,250:500])*3)
# #   #   arrows(wtstats[i,1,j-1],wtstats[i,2,j-1],wtstats[i,1,j],wtstats[i,2,j],lwd=4,
# #   #          col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),length=0.05)
# #   # }
# # }
# # frame()
# # legend("center",horiz=F,legend=strategies,pch=16,col=cols,bty="n",cex=2,ncol=2)
# 
# wtstats<-array(NA,dim=c(length(strategies),3,500))
# 
# for(i in 1:length(strategies))
# {
#   for(j in 1:500)
#   {
#     wtstats[i,1,j]<-wtd.mean(sp2trait[,j,i],weights=sp2cover[,j,i],normwt=T)
#     wtstats[i,2,j]<-wtd.var(sp2trait[,j,i],weights=sp2cover[,j,i],normwt=T)
#     wtstats[i,3,j]<-mean(sp2cover[,j,i])
#   }
# }
# par(mar=c(5,5,2,2)+.1)
# plot(1,1,type="n",xlim=c(26,36),ylim=c(0,20),bty="l",ylab="Weighted Variance Trait",
#      xlab="Weighted Mean Trait")
# 
# for(i in c(8,1,2,3,4,5,6,7))
# {
#   tempcols<-colorRampPalette(c("grey",cols[i]))(500)
#   afs<-wtstats[i,3,]/max(wtstats[,3,])
#   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16,cex=seq(1,500,by=1)/125)
#   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16)
#   for(j in 2:500)
#   {
#     points(wtstats[i,1,j],wtstats[i,2,j],
#            col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),
#            pch=16,cex=j/125)
#     arrows(wtstats[i,1,j-1],wtstats[i,2,j-1],wtstats[i,1,j],wtstats[i,2,j],lwd=4,
#            col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),length=0.05)
#   }
# }
# legend("topright",horiz=F,legend=strategies,pch=16,col=cols,bty="n",cex=1)
# 
# 
# for(i in c(8,1,2,3,4,5,6,7))
# {
#   print(i)
#   if(i==8) plot(density(diff(wtstats[i,2,])),col=adjustcolor(cols[i],alpha.f=.8),lwd=2,
#                 bty="n",main="")
#   else lines(density(diff(wtstats[i,2,])),col=adjustcolor(cols[i],alpha.f=.8),lwd=2)
#   
# }
# 
# install.packages("SciViews")
# library(SciViews)
# 
# xrange<-ceiling(max(abs((wtstats[,1,500]-wtstats[,1,1]))))
# yrange<-ceiling(max(abs((wtstats[,2,500]-wtstats[,2,1]))))
# 
# vecx<-(wtstats[,1,500]-wtstats[,1,1])/ceiling(max(abs((wtstats[,1,500]-wtstats[,1,1]))))
# vecy<-(wtstats[,2,500]-wtstats[,2,1])/ceiling(max(abs((wtstats[,2,500]-wtstats[,2,1]))))
# 
# 
# vectorplot(vecx,vecy,col=cols,lwd=2,ylab=expression(paste(Delta, " Trait ", Variance[w],sep=" ")),
#            xlab=expression(paste(Delta, " Trait ", Mean[w],sep=" ")),axes=F)
# axis(side=1,at=c(-1,-.5,0,.5,1),labels=c(-1*xrange,-.5*xrange,0,.5*xrange,xrange))
# axis(side=2,at=seq(-1,1,length.out=5),labels=c(-1*yrange,-.5*yrange,0,.5*yrange,yrange))
# legend("top",col=cols,legend=strategies,ncol=3,pch=16)
# # 
# par(mar=c(4,4,0,0))
# for(i in 1:7)
# {
#   plot(1,1,type="n",xlim=c(1,500),ylim=c(-.4,.4),bty="l",ylab=expression(paste(Delta, "Weighted Variance Trait")),xlab="Time")
#   aa<-diff(wtstats[i,2,])
#   tempcols<-colorRampPalette(c("grey",cols[i]))(500)
#   #afs<-wtstats[i,3,]/max(wtstats[,3,])
#   lines(seq(1,499),aa,col=tempcols[500],lwd=1)
#   polypoints<-rep(0,499*2)
# 
#   for(q in 1:(499))
#   {
#     polypoints[q]<-max(aa[q],0)
#   }
# 
#   for(q in (2*499):500)
#   {
#     print(q)
#     polypoints[q]<-min(aa[q-499],0)
#   }
#   bb<-cbind(c(seq(1,499,by=1),seq(499,1,by=-1)),polypoints)
#   polygon(bb, col=adjustcolor(tempcols[500],alpha.f=.8),border=F)
#   
#   abline(h=0,lwd=.5)
#   abline(h=mean(diff(wtstats[i,2,])),lty=2,lwd=2)
#   #print(mean(diff(wtstats[i,2,])))
#   abline(h=quantile(diff(wtstats[i,2,]),c(.025,.975))[c(1,2)])
#   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=afs),pch=16,cex=seq(1,500,by=1)/125)
#   # points(wtstats[i,1,],wtstats[i,2,],col=adjustcolor(tempcols,alpha.f=.5),pch=16)
#   # for(j in 2:500)
#   # {
#   #   
#   #   # points(wtstats[i,1,j],wtstats[i,2,j],col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),pch=16,cex=wtstats[i,3,j]/max(wtstats[,3,250:500])*3)
#   #   arrows(wtstats[i,1,j-1],wtstats[i,2,j-1],wtstats[i,1,j],wtstats[i,2,j],lwd=4,
#   #          col=adjustcolor(tempcols[j],alpha.f=wtstats[i,3,j]/max(wtstats[,3,])),length=0.05)
#   # }
# }
# frame()
# legend("center",horiz=F,legend=strategies,pch=16,col=cols,bty="n",cex=2,ncol=2)
# 

hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?
stratlabs<-c("Hot","Cold","Hot and Cold","High Cover","Low Cover","Gradient","Portfolio","Random")
windows() # opens graphics window
colramp<-colorRampPalette(c("blue","yellow","red"))(length(strategies))[c(2,1,seq(3,length(strategies),by=1))]
cols<-brewer.pal(9,"Paired")[-8]
#pdf("SingleSimulationExampleTimeSeries.pdf",onefile=T,height=7.6,width=10)
# layout(matrix(c(1,
#                 2,
#                 3,
#                 4),nrow=4,ncol=1))  # Layout panels according to nsp


pdf("SingleSimulationCoralCoverandTraitTrace6_diffs.pdf",height=6,width=6)
layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(1,1.25))
par(mar=c(1,5,1,1))
j<-1 # start counter

plot(seq(1,maxtime1),rep(0,maxtime1),ylim=c(-.1,.1),xlim=c(1500,maxtime1+500),
     col="black",type="l",lwd=2,yaxs="i",ylab="Coral Cover Diff. from Avg.",xaxt="n",xlab="",bty="l") # Plot the total hard coral density across the reef through time
# lines(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(size),])/size),ylim=c(0,1),xlim=c(0,2000),
#      col="darkorange",type="l",lwd=2) # Plot the total hard coral density across the reef through time
# 

for(i in 1:length(strategies))
{
  
  lines(seq(1,runtime)+maxtime1,(colSums(out2[((1-1)*size+1):(max(hard.corals)*size),,i])/size)-apply((colSums(out2[((1-1)*size+1):(max(hard.corals)*size),,])/size),MARGIN=1,FUN=mean),col=cols[i],lwd=1.2) # Plot the total hard coral density across the reef through time
  # lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(size),,i])/size),col=cols(7)[i],lwd=2,lty=2) # Plot the total hard coral density across the reef through time
  
}
legend("top",horiz=F,ncol=3,legend=stratlabs,lty=1,lwd=2,col=cols,bty="n",cex=.8)
legend("topleft","a",bty="n",cex=.8)
# dev.off()
# 
# pdf("SingleSimulationTraitTraces6.pdf",height=8,11)
par(mar=c(5,5,1,1))
plot(seq(1,maxtime1),rep(0,maxtime1),ylim=c(-1.5,1.0),xlim=c(1500,maxtime1+maxtime2),
     type="n",lwd=1.2,col="grey",lty=2,ylab="Optimal Temp. Diff. from Avg",xlab="Simulation Year",bty="l")
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(3*size),]
for(j in 1:maxtime1)
{
  
  
  
  weight.temps[,j]<-out$ts[(3*size+1):(5*size),j]*(out$ts[(0*size+1):(2*size),j]/(colSums(out$ts[(0*size+1):(2*size),])[j]))
  
  
}

lines(seq(1,maxtime1),rep(0,maxtime1),col="black",lty=1,lwd=1.2)

#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2)
for(i in 1:8)
{
  weight.tempsall<-out2[(1*size+1):(3*size),,1]
  for(j in 1:500)
  {
    
    #print((out2[(0*size+1):(2*size),j,i]/(colSums(out2[(0*size+1):(2*size),,i])[j])))
    #print(dim(mean(apply((out2[(0*size+1):(2*size),,]),MARGIN=3,FUN=colSums))[j,]))
    weight.tempsall[,j]<-out2[(3*size+1):(5*size),j,i]*(out2[(0*size+1):(2*size),j,i]/(colSums(out2[(0*size+1):(2*size),,i])[j]))-
      mean(out2[(3*size+1):(5*size),j,]*(out2[(0*size+1):(2*size),j,]/(mean(apply((out2[(0*size+1):(2*size),,]),MARGIN=3,FUN=colSums)[j,]))))
     
    #print(dim(weight.tempsall))
    
  }
  
  #   lines(seq(1,runtime)+1500,colMeans(out2[(3*size+1):(4*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  # lines(seq(1,runtime)+1500,colMeans(out$ts[(4*size+1):(5*size),,i]),col="dodgerblue",lty=2)
  lines(seq(1,runtime)+1500,colSums(weight.tempsall),col=cols[i],lty=1,lwd=1.2)
  
}
#legend("top",lty=2,col="grey",legend="Average Temperature",bty="n",lwd=1.2,cex=0.8)
#lines(seq(1,runtime)+maxtime1,colMeans(out2[(6*size+1):(7*size),,i]),xlim=c(1,maxtime1+maxtime2),lty=2,lwd=1.2,col="grey")
legend("topleft","b",bty="n",cex=.8)
dev.off()





pdf("SingleSimulationCoralCoverandTraitTrace6.pdf",height=6,width=6)
layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(1,1.25))
par(mar=c(1,5,1,1))
j<-1 # start counter

plot(seq(1,maxtime1),(colSums(out$ts[((1-1)*size+1):(max(hard.corals)*size),])/size),ylim=c(0,.8),xlim=c(1500,maxtime1+500),
     col="black",type="l",lwd=2,yaxs="i",ylab="Total Coral Cover",xaxt="n",xlab="",bty="l") # Plot the total hard coral density across the reef through time
# lines(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(size),])/size),ylim=c(0,1),xlim=c(0,2000),
#      col="darkorange",type="l",lwd=2) # Plot the total hard coral density across the reef through time
# 

for(i in 1:length(strategies))
{
  
  lines(seq(1,runtime)+maxtime1,(colSums(out2[((1-1)*size+1):(max(hard.corals)*size),,i])/size),col=cols[i],lwd=1.2) # Plot the total hard coral density across the reef through time
  # lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(size),,i])/size),col=cols(7)[i],lwd=2,lty=2) # Plot the total hard coral density across the reef through time
  
}
legend("top",horiz=F,ncol=3,legend=stratlabs,lty=1,lwd=2,col=cols,bty="n",cex=.8)
legend("topleft","a",bty="n",cex=.8)
# dev.off()
# 
# pdf("SingleSimulationTraitTraces6.pdf",height=8,11)
par(mar=c(5,5,1,1))
plot(seq(1,maxtime1),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(26,32),xlim=c(1500,maxtime1+maxtime2),
     type="l",lwd=1.2,col="grey",lty=2,ylab="Optimal Temperature",xlab="Simulation Year",bty="l")
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(3*size),]
for(j in 1:maxtime1)
{
  
  
  
  weight.temps[,j]<-out$ts[(3*size+1):(5*size),j]*(out$ts[(0*size+1):(2*size),j]/(colSums(out$ts[(0*size+1):(2*size),])[j]))
  
  
}

lines(seq(1,maxtime1),colSums(weight.temps),col="black",lty=1,lwd=1.2)

#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2)
for(i in 1:8)
{
  weight.temps<-out2[(1*size+1):(3*size),,i]
  for(j in 1:500)
  {
    
    
    
    weight.temps[,j]<-out2[(3*size+1):(5*size),j,i]*(out2[(0*size+1):(2*size),j,i]/(colSums(out2[(0*size+1):(2*size),,i])[j]))
    
    
  }
  
  #   lines(seq(1,runtime)+1500,colMeans(out2[(3*size+1):(4*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  # lines(seq(1,runtime)+1500,colMeans(out$ts[(4*size+1):(5*size),,i]),col="dodgerblue",lty=2)
  lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols[i],lty=1,lwd=1.2)
  
}
legend("top",lty=2,col="grey",legend="Average Temperature",bty="n",lwd=1.2,cex=0.8)
lines(seq(1,runtime)+maxtime1,colMeans(out2[(6*size+1):(7*size),,i]),xlim=c(1,maxtime1+maxtime2),lty=2,lwd=1.2,col="grey")
legend("topleft","b",bty="n",cex=.8)
dev.off()


pdf("SingleSimulationTraitTraces7.pdf",height=8,11)
layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(1,1.1))
par(mar=c(1,5,1,1))
plot(seq(1,maxtime1),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1,maxtime1+maxtime2),
     type="l",lwd=1.5,col="grey",lty=2,ylab="Optimal Temperature",xaxt="n",xlab="")
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:maxtime1)
{
  
  
  
  weight.temps[,j]<-out$ts[(3*size+1):(4*size),j]*(out$ts[(0*size+1):(1*size),j]/(colSums(out$ts[(0*size+1):(1*size),])[j]))
  
  
}

lines(seq(1,maxtime1),colSums(weight.temps),col="black",lty=1,lwd=1.5)

#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2)
for(i in 1:8)
{
  weight.temps<-out2[(1*size+1):(2*size),,i]
  for(j in 1:500)
  {
    
    

      weight.temps[,j]<-out2[(3*size+1):(4*size),j,i]*(out2[(0*size+1):(1*size),j,i]/(colSums(out2[(0*size+1):(1*size),,i])[j]))

    
  }
  lines(seq(1,runtime)+maxtime1,colMeans(out2[(6*size+1):(7*size),,i]),xlim=c(1,maxtime1+maxtime2),lty=2,lwd=1.5,col="grey")
#   lines(seq(1,runtime)+1500,colMeans(out2[(3*size+1):(4*size),,i]),col=cols(7)[i],lty=2,lwd=2)
 # lines(seq(1,runtime)+1500,colMeans(out$ts[(4*size+1):(5*size),,i]),col="dodgerblue",lty=2)
 lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols[i],lty=1,lwd=2)
  
}
legend("topleft",legend="Competitive",bty="n")

par(mar=c(5,5,1,1))
plot(seq(1,maxtime1),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1,maxtime1+maxtime2),
     type="l",lwd=1.5,col="grey",lty=2,ylab="Optimal Temperature",xlab="Year")
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2)
#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:1500)
{
  
  
  
  weight.temps[,j]<-out$ts[(4*size+1):(5*size),j]*(out$ts[(1*size+1):(2*size),j]/(colSums(out$ts[(1*size+1):(2*size),])[j]))
  
  
}



lines(seq(1,maxtime1),colSums(weight.temps),col="black",lty=1,lwd=1.5)

for(i in 1:9)
{
  weight.temps<-out2[(4*size+1):(5*size),,i]
  for(j in 1:500)
  {

      weight.temps[,j]<-out2[(4*size+1):(5*size),j,i]*(out2[(1*size+1):(2*size),j,i]/(colSums(out2[(1*size+1):(2*size),,i])[j]))

  
  }
  lines(seq(1,runtime)+1500,colMeans(out2[(6*size+1):(7*size),,i]),ylim=c(25,30),xlim=c(1,2000),
        lwd=1.5,lty=2,col="grey")
  #lines(seq(1,runtime)+1500,colMeans(out$ts[(3*size+1):(4*size),,i]),col="darkorange",lty=2)
  #lines(seq(1,runtime)+1500,colMeans(out2[(4*size+1):(5*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  
  lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols[i],lty=1,lwd=1.5)
}
legend("topleft",legend="Stress-Tolerant",bty="n")
legend("top",horiz=F,ncol=3,legend=stratlabs,lty=1,lwd=2,col=cols,bty="n",cex=.8)
dev.off()


image(mpa,col=c("turquoise1","lightcoral"))


layout(matrix(c(1,
                2,
                3,
                4),nrow=4,ncol=1))  # Layout panels according to nsp
par(mar=c(3,3,0,0))
j<-1 # start counter

plot(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(size),])/size),ylim=c(0,1),xlim=c(1250,2000),
     col="black",type="l",lwd=2,yaxs="i") # Plot the total hard coral density across the reef through time
# lines(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(size),])/size),ylim=c(0,1),xlim=c(0,2000),
#      col="darkorange",type="l",lwd=2) # Plot the total hard coral density across the reef through time
# 

for(i in 1:9)
{
  
  lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(size),,i])/size),col=cols[i],lwd=2) # Plot the total hard coral density across the reef through time
  # lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(size),,i])/size),col=cols(7)[i],lwd=2,lty=2) # Plot the total hard coral density across the reef through time
  
}
legend("top",horiz=T,legend=strategies,lty=1,lwd=2,col=cols,bty="n",cex=1.5)


plot(seq(1,maxtime),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1250,2000),
     type="l",lwd=2,col="grey",lty=2)
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:1500)
{
  
  
  
  weight.temps[,j]<-out$ts[(3*size+1):(4*size),j]*(out$ts[(1*size+1):(2*size),j]/(colSums(out$ts[(1*size+1):(2*size),])[j]))
  
  
}

lines(seq(1,maxtime),colSums(weight.temps),col="black",lty=1,lwd=2)

#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2)
for(i in 1:9)
{
  weight.temps<-out2[(1*size+1):(2*size),,i]
  for(j in 1:500)
  {
    
    
    
    weight.temps[,j]<-out2[(3*size+1):(4*size),j,i]*(out2[(1*size+1):(2*size),j,i]/(colSums(out2[(1*size+1):(2*size),,i])[j]))
    
    
  }
  lines(seq(1,runtime)+1500,colMeans(out2[(6*size+1):(7*size),,i]),xlim=c(1250,2000),lty=2,lwd=2,col="grey")
  #   lines(seq(1,runtime)+1500,colMeans(out2[(3*size+1):(4*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  # lines(seq(1,runtime)+1500,colMeans(out$ts[(4*size+1):(5*size),,i]),col="dodgerblue",lty=2)
  lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols[i],lty=1,lwd=2)
  
}

plot(seq(1,maxtime),(colSums(out$ts[((1)*size+1):(2*size),])/size),ylim=c(0,1),xlim=c(1250,2000),
     col="black",type="l",lwd=2,yaxs="i") # Plot the total hard coral density across the reef through time
# lines(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(size),])/size),ylim=c(0,1),xlim=c(0,2000),
#      col="darkorange",type="l",lwd=2) # Plot the total hard coral density across the reef through time
# 

for(i in 1:9)
{
  
  lines(seq(1,runtime)+1500,(colSums(out2[((1)*size+1):(2*size),,i])/size),col=cols[i],lwd=2) # Plot the total hard coral density across the reef through time
  # lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(size),,i])/size),col=cols(7)[i],lwd=2,lty=2) # Plot the total hard coral density across the reef through time
  
}

plot(seq(1,maxtime),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1250,2000),
     type="l",lwd=2,col="grey",lty=2)
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2)
#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:1500)
{
  
  
  
  weight.temps[,j]<-out$ts[(4*size+1):(5*size),j]*(out$ts[(2*size+1):(3*size),j]/(colSums(out$ts[(2*size+1):(3*size),])[j]))
  
  
}



lines(seq(1,maxtime),colSums(weight.temps),col="black",lty=1,lwd=2)

for(i in 1:9)
{
  weight.temps<-out2[(4*size+1):(5*size),,i]
  for(j in 1:500)
  {
    
    weight.temps[,j]<-out2[(4*size+1):(5*size),j,i]*(out2[(2*size+1):(3*size),j,i]/(colSums(out2[(2*size+1):(3*size),,i])[j]))
    
    
  }
  lines(seq(1,runtime)+1500,colMeans(out2[(6*size+1):(7*size),,i]),ylim=c(25,30),xlim=c(1250,2000),
        lty=2,lwd=2,col="grey")
  #lines(seq(1,runtime)+1500,colMeans(out$ts[(3*size+1):(4*size),,i]),col="darkorange",lty=2)
  #lines(seq(1,runtime)+1500,colMeans(out2[(4*size+1):(5*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  
  lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols[i],lty=1,lwd=2)
}

dev.off()


# add log ratio plot



















time.s<-out$ts
post.burn.sp1<-time.s[1:size,501:1000]
vcov.mat<-cov(t(post.burn.sp1))
dim(vcov.mat)
rowSums(vcov.mat)

portfolio.calc(t(post.burn.sp1),sites=c(3,8))

portfolio.calc<-function(x,sites)
{
  covx<-cov(x[,sites])
  Ex<-colMeans(x[,sites])
  x.vec<-rep(1,length(sites))/length(sites)
  sig.p.x2<-t(x.vec)%*%covx%*%x.vec
  mu.vec<-t(x.vec)%*%Ex
  return(list("ExpectedReturn"=mu.vec,"PortfolioVariance"=sig.p.x2))
}
possible.portfolios<-combn(25,5)
portfolio.perf<-matrix(NA,ncol=ncol(possible.portfolios),nrow=3)

for(i in 1:ncol(possible.portfolios))
{
  tempout<-portfolio.calc(t(post.burn.sp1),sites=possible.portfolios[,i])
  portfolio.perf[1,i]<-tempout$ExpectedReturn
  portfolio.perf[2,i]<-tempout$PortfolioVariance
  portfolio.perf[3,i]<-tempout$ExpectedReturn/sqrt(tempout$PortfolioVariance)
}
quick.port<-order(rowSums(vcov.mat))[c(1,2,3,4,5)]
full.port<-possible.portfolios[,which(portfolio.perf[2,]==min(portfolio.perf[2,]))]
print(quick.port)
print(full.port)
#===================================================================================
# Create plots of state values across space and time. These images will be 
#   heat maps for each state in the model. The color scales are not equal 
#   across panels.
#===================================================================================
imageCols<-colorRampPalette(c("blue","red","gold"))         # Generates a function to make a color ramp

layout.switch<-switch(as.character(nsp),"1"=matrix(c(1,3,            # Determines the plot layout dependent on the number of species in model
                                     2,3),nrow=2,ncol=2,byrow=T),
                      "2"=matrix(c(1,2,
                                 3,4,
                                 5,5),nrow=3,ncol=2,byrow=T),
                      "3"=matrix(c(1,2,3,
                                 4,5,6,
                                 0,7,0),nrow=3,ncol=3,byrow=T),
                      "4"=matrix(c(1,2,3,4,
                                 5,6,7,8,
                                 0,9,9,0),nrow=3,ncol=4,byrow=T))


all.labels<-c(paste("Coral",seq(1,nsp-1),"Density By Reef Location",sep=" "),"Macroalgae Density By Reef Location",
              paste("Coral",seq(1,nsp-1),"Trait By Reef Location",sep=" "),"Macroalgae Trait By Reef Location","Temperature")   # Vector of plot labels dependent on number of species
out3<-cbind(out$ts,out2[,,3])
windows() # opens graphics window
layout(layout.switch)  # generates panel layout according to layout.switch above
for(i in 1:(2*nsp+1))   # loop through each state variable
{
  if(i<=nsp){
    # Plot species densities (z) across time (x) and space (y)
    image(x=seq(1,2000),y=seq(1:size),z=t(out3[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
          ylab="Location on Reef",main=all.labels[i],zlim=c(0,1))
  }
  else{
    # Plot species traits or local temperature (z) across time (x) and space (y)
    image(x=seq(1,2000),y=seq(1:size),z=t(out3[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
          ylab="Location on Reef",main=all.labels[i],zlim=c(20,40))
  }
}

#====================================================================================
# Create plots of distribution of species densities over time
#   - 1 plot for hard corals
#   - 1 plot for macroalgae
#====================================================================================
hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?

cols<-c("dodgerblue","darkorange","darkgreen","purple","red")  # List of colors for line plots

layout.switch.2<-switch(as.character(length(macro.algae)>0),   # Switch for plot panel layout dependent on nsp
                        "TRUE"=matrix(c(1,2),nrow=2,ncol=1),
                        "FALSE"=matrix(c(1),nrow=1,ncol=1))

windows()  # opens graphics window
layout(layout.switch.2)  # Layout panels according to nsp

j<-1 # start counter
while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
{
  if(j==1) plot(seq(1,maxtime),(colSums(out$ts[((j-1)*size+1):(j*size),])/size),ylim=c(0,1),type="l",
               yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[j])                      
  else{lines(seq(1,maxtime),(colSums(out$ts[((j-1)*size+1):(j*size),])/size),col=cols[j])}
  j<-j+1
}
lines(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(max(hard.corals)*size),])/size),col="black",lwd=2) # Plot the total hard coral density across the reef through time

if(length(macro.algae)>0){ # Are there any macroalgae in the simulation
  while(j<=max(macro.algae)) # Loop through macroalgae species and plot the average density across the reef through time
  {
    if(j==min(macro.algae)) plot(seq(1,maxtime),(colSums(out$ts[((j-1)*size+1):(j*size),])/size),ylim=c(0,1),type="l",
                                yaxs="i",bty="l",ylab="MacroAlgae Cover",xlab="Time",col="black",lwd=2)
    else{lines(seq(1,maxtime),(colSums(out$ts[((j-1)*size+1):(j*size),])/size),col=cols[j])}
    j<-j+1
    
  }
}

#====================================================================================
# Create plots of relative contribution of ecological reorganization and 
#     evolutionary adaptation
#====================================================================================
rel_contrib<-dZbardt(Nall=out$ts[1:(nsp*size),],Zall=out$ts[(nsp*size+1):(2*nsp*size),]) # calculate ecological and evolutionary contributions

windows() # opens graphics window
plot(times[-1],log(abs(rel_contrib$Ecology/rel_contrib$Evolution)),col="darkgreen",type="l",ylim=c(-10,10)) # plot ecological turnover contribution across time
abline(h=0)
#lines(times[-1],rel_contrib$Evolution,col="darkorange") # plot evolutionary adaptation across time

#=====================================================================================
# Calculate Shannon Diversity Across Time
#=====================================================================================
spp1<-out$ts[1:size,] # Space in rows, time in columns
spp2<-out$ts[(size+1):(2*size),]
spp3<-out$ts[(size*2+1):(3*size),]

s.div<-matrix(0,nrow=size,ncol=maxtime)
shan.av<-rep(0,maxtime)
shan.max<-shan.av
shan.min<-shan.av
for(i in 1:maxtime)
{
  for(j in 1:size)
  {
#    subdat<-c(spp1[j,i],spp2[j,i],spp3[j,i])
    subdat<-c(spp1[j,i],spp2[j,i]) # don't include macroalgae
    tot.dens<-sum(subdat)
    s.div[j,i]<-shannon.div(subdat/tot.dens)
  }
  shan.av[i]<-mean(s.div[,i])
  shan.max[i]<-max(s.div[,i])
  shan.min[i]<-min(s.div[,i])
}

windows()
plot(times[-1],shan.av,type="l",lwd=2,ylim=c(0,1),yaxs="i",
     bty="l",ylab="Average Shannon Diversity",xlab="Time")
lines(times[-1],shan.max,lwd=1,col="grey30")
lines(times[-1],shan.min,lwd=1,col="grey30")

#====================================================================================
# Create plots of distribution of species densities over time
#   - 1 plot for hard corals
#   - 1 plot for macroalgae
#====================================================================================
hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?

cols<-c("dodgerblue","darkorange","darkgreen","purple","red")  # List of colors for line plots


windows()  # opens graphics window
#layout(layout.switch.2)  # Layout panels according to nsp

j<-1 # start counter
while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
{
  if(j==1) plot(seq(1,maxtime),(colSums(out$ts[((j+2)*size+1):((j+3)*size),])/size),ylim=c(25,35),type="l",
                yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[j])                      
  else{lines(seq(1,maxtime),(colSums(out$ts[((j+2)*size+1):((j+3)*size),])/size),col=cols[j])}
  j<-j+1
}
