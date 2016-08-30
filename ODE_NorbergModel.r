source("Norberg_Functions.r")

nsp<-2                        # how many species in model?
size<-20                      # how many reefs in model?
times<-seq(0,50,by=1)         # Vector from 1 to the number of time steps
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
sptype<-c(rep(1,(nsp-1)),1)
mpa<-matrix(rep(1,size*nsp),nrow=nsp,ncol=size,byrow=T)
mpa[sptype==2,c((round(size/3)):(2*round(size/3)))]<-1.5

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
  
  V=c(.1,.1),
  D=c(0.1,.1),
  rmax=c(1,.9),
  alphas=matrix(1,nrow=nsp,ncol=nsp),
  m=c(.1,.1),
  w=c(5,5),
  
  annual.temp.change=.02,
  maxtemp=30,
  Nmin=10^-6,
  deltax=1
  
  )


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
    temps<-y[(length(y)-(size-1)):length(y)]  # Extract current temperatures

    dspp<-spps  # create vector for change in densities
    dtraits<-traits  # create vector for change in traits
    
    for(i in 1:nsp)  # Loop through each species
    { 
      dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],
                     TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,]) # Calculate change in densities for species i across reef
      dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],
                        TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,],Nmin=Nmin)  # Calculate change in traits for species i across reef
    }
    
    # Calculate change in temperatures across reef
    dtemps<-switch(temp.change,
                   const=rep(0,size),
                   linear=rep(annual.temp.change,size),
                   sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size)) 

    dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector
    return(list(dsp))
  })

  
}

# Test model code for 1 time step
coral_trait(t=times,y=allstate,parms=parms,size=size,nsp=nsp,temp.change="const")

# Run model and return a matrix of state values across time
out <- ode.1D(y = allstate, t = times, func = coral_trait, parms = parms,
              nspec = (nsp*2)+1, names = allnames, dimens=size,size=size,nsp=nsp,temp.change="linear")

#===================================================================================
# Create plots of state values across space and time. These images will be 
#   heat maps for each state in the model. The color scales are not equal 
#   across panels.
#===================================================================================
delx<-1     # cell size for heat map grid
Distance <- seq(from = 0.5, by = delx, length.out = size) # Creates grid for heat map

windows() # open graphics window
spp.to.plot<-2 # How many species to plot?
which.to.plot<-c(1:spp.to.plot,((nsp+1):(nsp+spp.to.plot)),(nsp*2+1)) # Which states are being plotted. Alternatively can be set manually to picj individual states.

all.labels<-c(paste("Coral",seq(1,nsp),"Density By Reef Location",sep=" "),paste("Coral",seq(1,nsp),"Trait By Reef Location",sep=" "),"Temperature")
        # Creates vector of labels for heat map panels

image(out, which=allnames[which.to.plot], grid = Distance,
      xlab = "time", ylab = "Distance on Reef, km",
      main = all.labels[which.to.plot],legend=T)
