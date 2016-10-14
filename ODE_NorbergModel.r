source("Norberg_Functions.r")
start.parms<-read.csv("Norberg_StartingParameters_test.csv",header=T)

set.seed(3)
nsp<-3                        # how many species in model?
size<-50                      # how many reefs in model?
maxtime<-100                  # How long to run simulation?
times<-seq(0,maxtime,by=1)    # Vector from 1 to the number of time steps
mid<-25                       # mean temperature across all reefs at start of simulation.
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
sppstate<-generate.state(size=size,nsp=nsp,dens=.3,random=T)
traitstate<-generate.traits(nsp=nsp,size=size,mid=mid,range=range,temps=temps,trait.scenario="perfect.adapt")
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
mpa<-rep(0,size)
mpa[c((round(size/3)):(2*round(size/3)))]<-1

#==============================================================================
# Define model parameters (# required)
#   - nsp: Number of Species (1)
#   - V: Genetic variance (nsp)
#   - D: Dispersal (nsp)
#   - rmax: Maximum growth rate (nsp)
#   - alphas: species interaction (competition) matrix (nsp x nsp matrix)
#   - m: mortality rate (nsp)
#   - w: temperature tolerance or plasticity (nsp)
#   - annual.temp.change: rate of temperature change (1)
#   - maxtemp: maximum temperature (sigmoid temp change scenario) (1)
#   - Nmin: minimum density (for computational happiness) (1)
#   - deltax: step size in space for derivitives (1)
#==============================================================================

parms<-list(
  nsp=nsp,
  
  V=c(t(start.parms[start.parms[,1]=="V",((1:nsp)+1)])),
  D=c(t(start.parms[start.parms[,1]=="D",((1:nsp)+1)])),
  rmax=c(t(start.parms[start.parms[,1]=="Rmax",((1:nsp)+1)])),
  #alphas=matrix(1,nrow=nsp,ncol=nsp),
  #alphas=diag(1,nrow=nsp,ncol=nsp),
  alphas=matrix(c(1,1.1,1,
                  1,1,1,
                  .99,.99,1),nrow=nsp,ncol=nsp,byrow=T),
  m=c(t(start.parms[start.parms[,1]=="M.normal",((1:nsp)+1)])),

  w=c(t(start.parms[start.parms[,1]=="w",((1:nsp)+1)])),
  pcatastrophe=0.02,
  
  annual.temp.change=.02,
  maxtemp=30,
  Nmin=10^-6,
  deltax=1
  
)
# parms<-list(
#   nsp=nsp,
#   
#   V=c(.1,.1),
#   D=c(0.1,.1),
#   rmax=c(1,.9),
#   alphas=matrix(1,nrow=nsp,ncol=nsp),
#   m=c(.1,.1),
#   w=c(5,5),
#   
#   annual.temp.change=.02,
#   maxtemp=30,
#   Nmin=10^-6,
#   deltax=1
#   
#   )


#===================================================================================
# Function to calculate partial derivitives of population growth, trait change, and 
#     temperature change at a given time step. Returns a vector of value changes to
#     "allstate", defined above
#
# Parameters:
#     t: time step
#     y: Vector of state values (input "allstate" vector here)
#     parms: List of parameters to be used in the model
#     size: Size of reef
#     nsp: Number of species
#     names: vector of names for y
#     temp.change: Scenario of temperature change to be modeled. "const" has
#                   constant temperatures through time. "linear" has linear
#                   temperature change. "sigmoid" has logistic temperature change up
#                   to a maxtemp value.
#====================================================================================

coral_trait<-function(t,y, parms,size,nsp,names,temp.change=c("const","linear","sigmoid")){
  # Read in parameter values to use  
  with(as.list(parms), {
    
    spps<-matrix(y[1:(size*nsp)],nrow=nsp,ncol=size,byrow=T)  # Extract species density state values
    traits<-matrix(y[(size*nsp+1):(length(y)-size)],nrow=nsp,ncol=size,byrow=T) # Extract species trait state values
    temperatures<-y[(length(y)-(size-1)):length(y)]  # Extract current temperatures
    
    dspp<-spps  # create vector for change in densities
    dtraits<-traits  # create vector for change in traits
    
    for(i in 1:nsp)  # Loop through each species
    { 
      dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],
                     TC=temperatures,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],
                     growth="Huey") # Calculate change in densities for species i across reef
      dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],
                        TC=temperatures,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                        mortality.model="tempvary",spp=species[i],growth="Huey")  # Calculate change in traits for species i across reef
    }
    
    # Calculate change in temperatures across reef
    dtemps<-switch(temp.change,
                   const=rep(0,size),
                   linear=rep(annual.temp.change,size),
                   sigmoid=rep((annual.temp.change*mean(temperatures)*(1-(mean(temperatures)/maxtemp))),size)) 
    
    dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector
    return(list(dsp))
  })
  
  
}

# Test model code for 1 time step
coral_trait(t=times,y=allstate,parms=parms,size=size,nsp=nsp,temp.change="sigmoid")

# Run model and return a matrix of state values across time
out <- ode.1D(y = allstate, t = times, func = coral_trait, parms = parms,
              nspec = (nsp*2)+1, names = allnames, dimens=size,size=size,nsp=nsp,temp.change="sigmoid")

#===================================================================================
# Create plots of state values across space and time. These images will be 
#   heat maps for each state in the model. The color scales are not equal 
#   across panels.
#===================================================================================
delx<-1     # cell size for heat map grid
Distance <- seq(from = 0.5, by = delx, length.out = size) # Creates grid for heat map

windows() # open graphics window
spp.to.plot<-3 # How many species to plot?
which.to.plot<-c(1:spp.to.plot,((nsp+1):(nsp+spp.to.plot)),(nsp*2+1)) # Which states are being plotted. Alternatively can be set manually to picj individual states.

all.labels<-c(paste("Coral",seq(1,nsp),"Density By Reef Location",sep=" "),paste("Coral",seq(1,nsp),"Trait By Reef Location",sep=" "),"Temperature")
# Creates vector of labels for heat map panels

image(out, which=allnames[which.to.plot], grid = Distance,
      xlab = "time", ylab = "Distance on Reef, km",
      main = all.labels[which.to.plot],legend=T)


rel_contrib<-dZbardt(Nall=t(out[,2:(nsp*size+1)]),Zall=t(out[,(nsp*size+2):(2*nsp*size+1)])) # calculate ecological and evolutionary contributions

windows() # opens graphics window
plot(times,rel_contrib$Ecology,col="darkgreen",type="l",ylim=c(0,.1)) # plot ecological turnover contribution across time
lines(times,rel_contrib$Evolution,col="darkorange") # plot ecological turnover contribution across time

abline(h=0)

#====================================================================================
# Create plots of distribution of species densities over time
#   - 1 plot for hard corals
#   - 1 plot for macroalgae
#====================================================================================
hard.corals<-which(species=="C") # Which of the species are hard corals?
macro.algae<-which(species=="MA") # Which of the species are macroalgae?

cols<-c("dodgerblue","darkorange","darkgreen","purple","red")  # List of colors for line plots

layout.switch.2<-switch(as.character(length(macro.algae)>0),   # Switch for plot panel layout dependent on nsp
                        "TRUE"=matrix(c(1,2),nrow=2,ncol=1),
                        "FALSE"=matrix(c(1),nrow=1,ncol=1))

windows()  # opens graphics window
layout(layout.switch.2)  # Layout panels according to nsp

j<-1 # start counter
while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
{
  if(j==1) plot(seq(1,(maxtime+1)),(rowSums(out[,(((j-1)*size+1)+1):((j*size)+1)]))/size,ylim=c(0,1),type="l",
                yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time",col=cols[j])                      
  else{lines(seq(1,(maxtime+1)),(rowSums(out[,(((j-1)*size+1)+1):((j*size)+1)])/size),col=cols[j])}
  j<-j+1
}
lines(seq(1,(maxtime+1)),(rowSums(out[,(((1-1)*size+1)+1):((max(hard.corals)*size)+1)])/size),col="black",lwd=2) # Plot the total hard coral density across the reef through time

if(length(macro.algae)>0){ # Are there any macroalgae in the simulation
  while(j<=max(macro.algae)) # Loop through macroalgae species and plot the average density across the reef through time
  {
    if(j==min(macro.algae)) plot(seq(1,(maxtime+1)),(rowSums(out[,(((j-1)*size+1)+1):((j*size)+1)])/size),ylim=c(0,1),type="l",
                                 yaxs="i",bty="l",ylab="MacroAlgae Cover",xlab="Time",col="black",lwd=2)
    else{lines(seq(1,(maxtime+1)),(rowSums(out[,(((j-1)*size+1)+1):((j*size)+1)])/size),col=cols[j])}
    j<-j+1
    
  }
}

