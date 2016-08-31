source("Norberg_Functions.r")

nsp<-3                        # how many species in model?
size<-20                      # how many reefs in model?
maxtime<-50                   # how many time steps?
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
sppstate<-generate.state(size=size,nsp=nsp,dens=0.1)
traitstate<-generate.traits(nsp=nsp,size=size,mid=mid,range=range,trait.scenario="u.const")
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
sptype<-c(rep(1,(nsp-1)),2)
mpa<-matrix(rep(1,size*nsp),nrow=nsp,ncol=size,byrow=T)
mpa[sptype==2,c((round(size/3)):(2*round(size/3)))]<-1.2

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
  
  V=c(.1,.1,.1),
  D=c(0.1,.1,.1),
  rmax=c(1,.9,.8),
  alphas=matrix(1,nrow=nsp,ncol=nsp),
  m.normal=c(.1,.1,.1),
  m.catastrophe=c(.7,.15,.05),
  w=c(5,5,5),
  pcatastrophe=0.05,
  
  annual.temp.change=.02,
  maxtemp=30,
  Nmin=10^-6,
  deltax=1
  
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
      temps<-time.series[(length(y)-(size-1)):length(y),k-1]  # Extract temperature state values
      
      pcat<-rbinom(1,1,pcatastrophe)           # Binomial draw for annual catastrophe occurance
      m<-(1-pcat)*m.normal+pcat*m.catastrophe  # set annual mortality rate dependent on presence of catastrophe
      
      dspp<-spps   # create storage for change in density state values
      dtraits<-traits # create storage for change in trait state values
      
      for(i in 1:nsp)
      { 
        dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,])
        dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,],Nmin=Nmin)
      }
      
      # Calculate change in temperatures across reef
      dtemps<-switch(temp.change,
                     const=rep(0,size),
                     linear=rep(annual.temp.change,size),
                     sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))+stoch.temp*rnorm(size,0,.1)
      
      dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector
      

      dspout[,k]<-dsp  # Store change in state value
      time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
      
    }
    return(list("ts"=time.series,"dsp"=dspout))
  })


  
}

out<-coral_trait_stoch(t=maxtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                       temp.change="const")

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
                                 9,9,10,11),nrow=3,ncol=4,byrow=T))


all.labels<-c(paste("Coral",seq(1,nsp),"Density By Reef Location",sep=" "),
              paste("Coral",seq(1,nsp),"Trait By Reef Location",sep=" "),"Temperature")   # Vector of plot labels dependent on number of species

windows() # opens graphics window
layout(layout.switch)  # generates panel layout according to layout.switch above
for(i in 1:(2*nsp+1))   # loop through each state variable
{
  if(i<=nsp){
    # Plot species densities (z) across time (x) and space (y)
    image(x=times,y=seq(1:size),z=t(out$ts[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
          ylab="Location on Reef",main=all.labels[i],zlim=c(0,.9))
  }
  else{
    # Plot species traits or local temperature (z) across time (x) and space (y)
    image(x=times,y=seq(1:size),z=t(out$ts[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
          ylab="Location on Reef",main=all.labels[i],zlim=c(22,28))
  }
}

#====================================================================================
# Create plots of distribution of species densities over time
#   - 1 plot for hard corals
#   - 1 plot for macroalgae
#
# NOT CURRENTLY AUTOMATED FOR CHANGES IN NSP!!!
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

