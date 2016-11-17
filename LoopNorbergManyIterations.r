# Loop stochastic Norberg Model across parameter and management scenarios
#=========================================================================
# Read in data files and source functions
source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1) # Matrix of parameter values to be filled from "scenarios"
scenarios<-read.csv("Norberg_StartingParameters_scenariolist.csv",header=T) # Each row contains parameter values for one scenario. Column headers identify parameter.

# Set initialization parameters
nsp<-3                        # how many species in model?
size<-60                      # how many reefs in model?
burnin<-1000                  # How long to run burn-in period?
runtime<-500                 # How long to run management/environmental change period?
mid<-28                       # mean temperature across all reefs at start of simulation.
range<-5                      # range of temperatures across reefs at start of simulation
iterations<-1                 # Number of stochastic iterations
nstrat<-length(unique(scenarios[,ncol(scenarios)])) # Number of management strategies
monitor.yrs<-c(50,100,500,1000)

# Create Multivariate normal covariance matrix for temperature anomalies
mdim<-size
lindec<-exp(seq(0,-5,length=mdim))# Decrease the second value to reduce the range of correlation
ma<-matrix(0,nrow=mdim,ncol=mdim)
ma[,1]<-lindec
for(i in 2:mdim){
  ma[-c(1:(i-1)),i]<-lindec[1:(mdim-(i-1))]
}
ma<-ma+t(ma)-diag(nrow(ma))
is.positive.definite(ma) # Is the matrix positive definite? Needs to be.

sds<-rep(.24,size) # SD of intersite temperature variation
b<-sds%*%t(sds)
spatialtemp<-b*ma # Scales MVN matrix to 0.2 * temp.stoch when temp.stoch=1.2


# Run model iterations in parallel across all cores on computer
cl<-makeCluster(detectCores())
registerDoParallel(cl)

tsout<-foreach(qqq=1:nrow(scenarios),.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{
#tsout<-foreach(qqq=42,.packages=c("deSolve","numDeriv","mvtnorm")) %do%{
  start.parms<-matrix(scenarios[qqq,][-c(1,ncol(scenarios))],nrow=5,ncol=3,byrow=F) # fill start.parms matrix from scenarios file
  colnames(start.parms)<-colnames(start.parm.mat) # identify columns
  rownames(start.parms)<-rownames(start.parm.mat) # identify rows
  
  strategies<-c("hot","cold","highcoral","lowcoral","portfolio","random","none") # Name different management strategies
  a<-runMod(q=qqq,strategy=strategies[scenarios[qqq,ncol(scenarios)]],
            iters=iterations,burnin=burnin,runtime=runtime,spatialtemp=spatialtemp,
            monitor.yrs=monitor.yrs) # Run the model for given scenario
  a
}

registerDoSEQ() # Close out cluster/parallel operations

#save.image(ADD LOCATION)


#=================================================================
# Summary Figures
#=================================================================
species<-c("C1","C2","MA") # Species Names
sptype<-c(1,1,2) # Type of species (Coral = 1, Macroalgae = 2)
hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?

cols<-c("dodgerblue","darkorange","darkgreen")  # List of colors for line plots

layout.switch.2<-switch(as.character(length(tsout)/nstrat),   # Switch for plot panel layout dependent on nsp
                        "7"=matrix(c(1,2,3,
                                     4,5,6,
                                     0,7,0,
                                     8,8,8),nrow=4,ncol=3,byrow=T),
                        "8"=matrix(c(1,2,3,
                                     4,5,6,
                                     7,0,8,
                                     9,9,9),nrow=4,ncol=3,byrow=T),
                        "9"=matrix(c(1,2,3,
                                     4,5,6,
                                     7,8,9,
                                     10,10,10),nrow=4,ncol=3,byrow=T))

windows()  # opens graphics window
layout(layout.switch.2,heights=c(1,1,1,.5))  # Layout panels according to nsp
for(i in 1:9)
{

  panel.use<-seq(i,i+(9*6),by=9) # Select simulation results that have the same ecological parameters (but different management strategies)
  #print(panel.use)
  for(z in 1:length(panel.use))
  {
    #print(panel.use[z])
    j<-1 # start counter
    while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
    {
      x.length<-length((colSums(tsout[[panel.use[z]]]$ExampleTimeSeries[((j-1)*size+1):(j*size),])/size)) # How long is the time series?
      if(j==1&z==1) {plot(seq(1,x.length),(colSums(tsout[[panel.use[z]]]$ExampleTimeSeries[((j-1)*size+1):(j*size),])/size),ylim=c(0,1),type="l",
                         yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[j])       } # Start new plotting panel for Coral 1           
      else{lines(seq(1,x.length),(colSums(tsout[[panel.use[z]]]$ExampleTimeSeries[((j-1)*size+1):(j*size),])/size),col=cols[j])} # Add Coral 2 to existing panel
      j<-j+1 # increase counter value
    }
    lines(seq(1,x.length),(colSums(tsout[[panel.use[z]]]$ExampleTimeSeries[((1-1)*size+1):(max(hard.corals)*size),])/size),col="black",lwd=2) # Plot the total hard coral density across the reef through time
    
    if(length(macro.algae)>0){ # Are there any macroalgae in the simulation
      while(j<=max(macro.algae)) # Loop through macroalgae species and plot the average density across the reef through time (on same panel as coral species)
      {
        lines(seq(1,x.length),(colSums(tsout[[panel.use[z]]]$ExampleTimeSeries[((j-1)*size+1):(j*size),])/size),col=cols[j],lwd=2)
        j<-j+1
        
      }
    }
  }
  mtext(side=3,paste("C.V",scenarios[i,2],"C.D",scenarios[i,3],"T.V",scenarios[i,7],"T.D",scenarios[i,8])) # label each panel with ecological parameters
}
frame() # Create blank panel for legend
legend("center",horiz=T,bty="n",lty=1,col=c(cols,"black"),legend=c("Competitive","Stress-Tolerant","Macroalgae",
                                                   "Hard Corals (Sum)"), cex=2) # Add legend

# Plot hard coral cover by management strategy for each parameter set
windows()  # opens new graphics window
layout(layout.switch.2,heights=c(1,1,1,.5))  # Layout panels according to nsp
strategies<-c("hot","cold","highcoral","lowcoral","portfolio","random","none") # identify strategies. Order corresponds to last column in scenarios
cols<-c("dodgerblue","dodgerblue","darkorange","darkorange","darkgreen","darkgreen","grey") 
linetypes<-c(rep(c(1,2),3),1)
for(i in 1:9)
{
  panel.use<-seq(i,i+(9*6),by=9)
  for(z in 1:length(panel.use))
  {
    x.length<-length((colSums(tsout[[panel.use[z]]][((1-1)*size+1):(1*size),,1])/size))
      if(z==1) plot(seq(1,x.length),(colSums(tsout[[panel.use[z]]][((1-1)*size+1):(max(hard.corals)*size),,1])/size),ylim=c(0,1),type="l",
                         yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[z],lty=linetypes[z],lwd=2)                      
      else{lines(seq(1,x.length),(colSums(tsout[[panel.use[z]]][((1-1)*size+1):(max(hard.corals)*size),,1])/size),col=cols[z],lwd=2,lty=linetypes[z])} # Plot the total hard coral density across the reef through time

  }
 mtext(side=3,paste("C.V",scenarios[i,2],"C.D",scenarios[i,3],"T.V",scenarios[i,7],"T.D",scenarios[i,8]))
        
}
frame()
legend("center",horiz=T,lty=linetypes,col=cols,legend=strategies,cex=2,lwd=2,bty="n")
  
  


coralsums<-colSums(tsout[[4]][((2-1)*size+1):(2*size),])/size

quantile(coralsums,0.975)

fitdist(coralsums,"beta")
fitlogit<-twCoefLogitnorm(median=median(coralsums),quant=quantile(coralsums,0.975),perc=0.975)

hist(rlogitnorm(10000,mu=fitlogit[1],sigma=fitlogit[2]))
hist(sample(coralsums,size=10000,replace=T),add=T,col="grey")





runMod<-function(q,strategy,iters,burnin,runtime,spatialtemp,monitor.yrs)
{ 
  nsp<-3                        # how many species in model?
  size<-60                      # how many reefs in model?
  maxtime<-runtime                  # how many time steps?
  times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
  mid<-28                       # mean temperature across all reefs at start of simulation.
  range<-5    # range of temperatures across reefs at start of simulation
  #print("Yo")
  stoch.time.series<-matrix(data=NA,nrow=(2*nsp*size+size),ncol=maxtime) # Create storage matrix for time series of state values at each location
  output.list<-list("ExampleTimeSeries"=stoch.time.series,"C1"=matrix(NA,nrow=iters,ncol=(length(monitor.yrs)*2)),"C2"=matrix(NA,nrow=iters,ncol=(length(monitor.yrs)*2)),
                    "Call"=matrix(NA,nrow=iters,ncol=(length(monitor.yrs)*2)),"MA"=matrix(NA,nrow=iters,ncol=(length(monitor.yrs)*2)),"Shannon"=matrix(NA,nrow=iters,ncol=(length(monitor.yrs)*2)))
  #  print("A")
  for(seednum in 1:iters)
  {
    #   set.seed(seeds[q])
    set.seed(seednum)
    #set.seed(1)
    nsp<-3                        # how many species in model?
    size<-60                      # how many reefs in model?
    maxtime1<-burnin              # how many time steps?
    times<-seq(0,maxtime1,by=1)   # Vector from 1 to the number of time steps
    mid<-28                       # mean temperature across all reefs at start of simulation.
    range<-5                      # range of temperatures across reefs at start of simulation
    #print("B")
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
    print("C")
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
    
    mpa<-setMPA(temps=temps,Nall=matrix(sppstate,nrow=nsp,ncol=size,byrow=T),sptype=sptype,size=size,amount=0.2,
                strategy="none")


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
      deltax=1,
      spatialtemp=spatialtemp
      
    )
    
    #print("D")
    out1<-coral_trait_stoch(t=maxtime1,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                            temp.change="const",stoch.temp=T)
    print("BurnCheck")
    #print(out$ts[1,maxtime])
    maxtime<-runtime                   # how many time steps?
    times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps

    temps<-generate.temps(size=size,mid=mid,range=range,temp.scenario="linear") 
    sppstate<-out1$ts[(1:(nsp*size)),maxtime1]
    traitstate<-out1$ts[(nsp*size+1):(2*nsp*size),maxtime1]
    allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)
    
    allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")


    
    mpa<-setMPA(temps=temps,Nall=matrix(sppstate,nrow=nsp,ncol=size,byrow=T),sptype=sptype,size=size,amount=0.2,
                strategy=strategy,priordata=out1$ts)
    print("MPACheck")
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
      # alphas=matrix(1,nrow=nsp,ncol=nsp),
      # alphas=diag(1,nrow=nsp,ncol=nsp),
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
      deltax=1,
      spatialtemp=spatialtemp
    )
    
    #print("E")
    out<-coral_trait_stoch(t=maxtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                           temp.change="sigmoid",stoch.temp=T)
    if(seednum==1) output.list$ExampleTimeSeries[,]<-out$ts
       
    c1means<-colMeans(out$ts[1:size,])[monitor.yrs]
    c2means<-colMeans(out$ts[(1+size):(2*size),])[monitor.yrs]
    callmeans<-c1means+c2means
    mameans<-colMeans(out$ts[(2*size+1):(3*size),])[monitor.yrs]
    
    c1sd<-apply(out$ts[1:size,],MARGIN=2,FUN=sd)[monitor.yrs]
    c2sd<-apply(out$ts[(1+size):(2*size),],MARGIN=2,FUN=sd)[monitor.yrs]
    callsd<-apply((out$ts[1:(size),]+out$ts[(size+1):(2*size),]),MARGIN=2,FUN=sd)[monitor.yrs]
    masd<-apply(out$ts[(2*size+1):(3*size),],MARGIN=2,FUN=sd)[monitor.yrs]
    
    output.list$C1[seednum,seq(1,2*length(monitor.yrs)-1,by=2)]<-c1means
    output.list$C2[seednum,seq(1,2*length(monitor.yrs)-1,by=2)]<-c2means
    output.list$Call[seednum,seq(1,2*length(monitor.yrs)-1,by=2)]<-callmeans
    output.list$MA[seednum,seq(1,2*length(monitor.yrs)-1,by=2)]<-mameans
    
    output.list$C1[seednum,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-c1sd
    output.list$C2[seednum,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-c2sd
    output.list$Call[seednum,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-callsd
    output.list$MA[seednum,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-masd
  }
  return(output.list)
}

#===================================================================================
# Function to calculate partial derivitives of population growth, trait change, and 
#     temperature change at a given time step. Appends the new state values to 
#     a matrix storing a time series of state values. Additionally creates a time 
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
    anoms<-matrix(NA,nrow=size,ncol=(t-1))
    dspout<-time.series # Create storage matrix for time series of change in state values at each location
    time.series[,1]<-y  # set initial state
    
    for(k in 2:t)
    {
      print(k) # Print counter
      
      
      spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
      traits<-matrix(time.series[(size*nsp+1):(length(y)-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
      anoms[,k-1]<-stoch.temp*rnorm(1,0,temp.stoch)+stoch.temp*rmvnorm(1,rep(0,size),spatialtemp)
      temps<-time.series[(length(y)-(size-1)):length(y),k-1]+anoms[,k-1] # Extract temperature state values
      #print("Anomalies Check")
#       print(temps)
#       print(anoms)
      dspp<-spps   # create storage for change in density state values
      dtraits<-traits # create storage for change in trait state values
      

      for(i in 1:nsp)
      { 
        #print("AA")
        dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                       mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
        #print("Test")
        dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                          mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                          mortality.model="tempvary",spp=species[i],growth="normal")
        #print("test2")
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
    return(list("ts"=time.series,"dsp"=dspout,"TemperatureAnomalies"=anoms))
  })
  
  
  
}
