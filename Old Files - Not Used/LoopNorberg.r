# Loop stochastic Norberg
#======================================

source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1)
scenarios<-read.csv("Norberg_StartingParameters_scenariolist.csv",header=T)

library(doSNOW)
library(doParallel)
library(foreach)
library(rbenchmark)
library(fitdistrplus)
library(logitnorm)
seeds<-seq(1,5,by=1)
nsp<-3                        # how many species in model?
size<-60                      # how many reefs in model?
maxtime<-1000                 # how many time steps?
times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
mid<-28                       # mean temperature across all reefs at start of simulation.
range<-5                      # range of temperatures across reefs at start of simulation

tsout<-matrix(NA,nrow=5,ncol=((nsp*size*2)+size))

cl<-makeCluster(detectCores())
registerDoParallel(cl)

tsout<-foreach(qqq=1:54,.packages=c("deSolve","numDeriv")) %dopar%{
  
  start.parms<-matrix(scenarios[qqq,][-c(1,ncol(scenarios))],nrow=5,ncol=3,byrow=F)
  print("A")
  strategies<-c("hot","cold","highcoral","lowcoral","random","none")
  print("B")
  colnames(start.parms)<-colnames(start.parm.mat)
  print("C")
  rownames(start.parms)<-rownames(start.parm.mat)
  print("D")
  a<-runMod(q=qqq,strategy=strategies[scenarios[qqq,ncol(scenarios)]])
  a
}

registerDoSEQ()

species<-c("C1","C2","MA")
sptype<-c(1,1,2)
hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?

cols<-c("dodgerblue","darkorange","darkgreen")  # List of colors for line plots

layout.switch.2<-switch(as.character(length(tsout)/6),   # Switch for plot panel layout dependent on nsp
                        "7"=matrix(c(1,2,3,
                                     4,5,6,
                                     0,7,0),nrow=3,ncol=3,byrow=T),
                        "8"=matrix(c(1,2,3,
                                     4,5,6,
                                     7,0,8),nrow=3,ncol=3,byrow=T),
                        "9"=matrix(c(1,2,3,
                                     4,5,6,
                                     7,8,9,
                                     10,10,10),nrow=4,ncol=3,byrow=T))

windows()  # opens graphics window
layout(layout.switch.2,heights=c(1,1,1,.5))  # Layout panels according to nsp
for(i in 1:9)
{

  panel.use<-seq(i,i+(9*5),by=9)
  print(panel.use)
  for(z in 1:length(panel.use))
  {
    print(panel.use[z])
    j<-1 # start counter
    while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
    {
      if(j==1&z==1) {plot(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((j-1)*size+1):(j*size),])/size),ylim=c(0,1),type="l",
                         yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[j])       }               
      else{lines(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((j-1)*size+1):(j*size),])/size),col=cols[j])}
      j<-j+1
    }
    lines(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((1-1)*size+1):(max(hard.corals)*size),])/size),col="black",lwd=2) # Plot the total hard coral density across the reef through time
    
    if(length(macro.algae)>0){ # Are there any macroalgae in the simulation
      while(j<=max(macro.algae)) # Loop through macroalgae species and plot the average density across the reef through time
      {
        lines(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((j-1)*size+1):(j*size),])/size),col=cols[j],lwd=2)
        j<-j+1
        
      }
    }
  }
  mtext(side=3,paste("C.V",scenarios[i,2],"C.D",scenarios[i,3],"T.V",scenarios[i,7],"T.D",scenarios[i,8]))
}
frame()
legend("center",horiz=T,bty="n",lty=1,col=c(cols,"black"),legend=c("Competitive","Stress-Tolerant","Macroalgae",
                                                   "Hard Corals (Sum)"), cex=2)


windows()  # opens graphics window
layout(layout.switch.2,heights=c(1,1,1,.5))  # Layout panels according to nsp
cols<-c("dodgerblue","dodgerblue","darkorange","darkorange","red","black")
linetypes<-rep(c(1,2),3)
for(i in 1:9)
{
  panel.use<-seq(i,i+(9*5),by=9)
  for(z in 1:length(panel.use))
  {
      
      if(z==1) plot(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((1-1)*size+1):(max(hard.corals)*size),])/size),ylim=c(0,2),type="l",
                         yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[z],lty=linetypes[z],lwd=2)                      
      else{lines(seq(1,maxtime),(colSums(tsout[[panel.use[z]]][((1-1)*size+1):(max(hard.corals)*size),])/size),col=cols[z],lwd=2,lty=linetypes[z])} # Plot the total hard coral density across the reef through time

  }
 mtext(side=3,paste("C.V",scenarios[i,2],"C.D",scenarios[i,3],"T.V",scenarios[i,7],"T.D",scenarios[i,8]))
        
}
frame()
legend("center",horiz=T,lty=linetypes[c(1,2,3,4,5,6)],col=cols[c(1,2,3,4,5,6)],legend=strategies[c(1,2,3,4,5,6)],cex=2,lwd=2,bty="n")
  
  


coralsums<-colSums(tsout[[4]][((2-1)*size+1):(2*size),])/size

quantile(coralsums,0.975)

fitdist(coralsums,"beta")
fitlogit<-twCoefLogitnorm(median=median(coralsums),quant=quantile(coralsums,0.975),perc=0.975)

hist(rlogitnorm(10000,mu=fitlogit[1],sigma=fitlogit[2]))
hist(sample(coralsums,size=10000,replace=T),add=T,col="grey")

runMod<-function(q,strategy)
{
  
#   set.seed(seeds[q])
  set.seed(2)
  nsp<-3                        # how many species in model?
  size<-60                      # how many reefs in model?
  maxtime<-1000                   # how many time steps?
  times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
  mid<-28                       # mean temperature across all reefs at start of simulation.
  range<-5                      # range of temperatures across reefs at start of simulation
  
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

  mpa<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=0.2,
              strategy=strategy)
  #print(mpa)
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
    mpa=mpa,
    species=species,
    sptype=sptype,
    V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
    D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
    rmax=as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
    #  alphas=matrix(1,nrow=nsp,ncol=nsp),
    #alphas=diag(1,nrow=nsp,ncol=nsp),
    alphas=matrix(c(1,1.1,1,
                    1,1,1,
                    .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
    m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
    
    w=c(t(as.numeric(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
    pcatastrophe=0.02,
    temp.stoch=1.2, #0.9 works great for constant temp
    annual.temp.change=.011,
    maxtemp=30,
    Nmin=10^-6,
    deltax=1
#     nsp=nsp,
#     mpa=mpa,
#     species=species,
#     sptype=sptype,
#     V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
#     D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
#     rmax=as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
#     #  alphas=matrix(1,nrow=nsp,ncol=nsp),
# #    alphas=diag(1,nrow=nsp,ncol=nsp),
#      
#      alphas=matrix(c(1,1.1,1,
#                      1,1,1,
#                      .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
#     m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
#     w=as.numeric(c(t(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
#     pcatastrophe=0.02,
#     temp.stoch=1.2, #0.9 works great for constant temp
#     annual.temp.change=.011,
#     maxtemp=30,
#     Nmin=10^-6,
#     deltax=1
    
  )

  
  out<-coral_trait_stoch(t=maxtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                         temp.change="sigmoid",stoch.temp=T)
  print(out$ts[,maxtime])
  return(out$ts)
}

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
                            stoch.temp=TRUE){
  
  with(as.list(parms), {  # Load parameter values
    time.series<-matrix(NA,nrow=(2*nsp*size+size),ncol=t) # Create storage matrix for time series of state values at each location
    dspout<-time.series # Create storage matrix for time series of change in state values at each location
    time.series[,1]<-y  # set initial state
    
    for(k in 2:t)
    {
      print(k) # Print counter
      
      
      spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
      traits<-matrix(time.series[(size*nsp+1):(length(y)-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
      temps<-time.series[(length(y)-(size-1)):length(y),k-1]+stoch.temp*rnorm(size,0,temp.stoch) # Extract temperature state values
      
      print("Yo")
      #       
      #       pcat<-rbinom(1,1,pcatastrophe)           # Binomial draw for annual catastrophe occurance
      #       m<-(1-pcat)*m.normal+pcat*m.catastrophe  # set annual mortality rate dependent on presence of catastrophe
      #       

      dspp<-spps   # create storage for change in density state values
      dtraits<-traits # create storage for change in trait state values
      
      print("Hey")
      print(rmax)
      print(mpa)
      print(species)
      for(i in 1:nsp)
      { 
        print("A")
        dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                       mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
        print("Test")
        dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                          mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                          mortality.model="tempvary",spp=species[i],growth="normal")
        print("test2")
      }
      
      # Calculate change in temperatures across reef
      dtemps<-switch(temp.change,
                     const=rep(0,size),
                     linear=rep(annual.temp.change,size),
                     sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))
      
      dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector
      
      
      dspout[,k]<-dsp  # Store change in state value
      time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
      
    }
    return(list("ts"=time.series,"dsp"=dspout))
  })
  
  
  
}
