source("Norberg_Functions.r")
library(mvtnorm)
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1)
scenarios<-read.csv("Norberg_StartingParameters_iterationlist_reservesize.csv",header=T)

start.parms<-matrix(scenarios[5,-c(1,17,18)],nrow=5,ncol=3,byrow=F)
colnames(start.parms)<-colnames(start.parm.mat)
rownames(start.parms)<-rownames(start.parm.mat)


cl<-makeCluster(detectCores()-1) # Set up parallel processing cluster
registerDoParallel(cl)
start.time<-Sys.time() # Save start time mark

#===============================================================================================
# Run the simulation model in parallel using foreach. This will run 100 iterations of the model
#===============================================================================================
mod.outall<-foreach(qqq=1:100,.packages=c("deSolve","numDeriv","mvtnorm","matrixcalc")) %dopar% {
  set.seed(qqq)                 # Set new seed value for random number generation
  nsp<-3                        # how many species in model?
  size<-60                      # how many reefs in model?
  maxtime1<-1000                # Length of burn in period
  burnin<-maxtime1              # Length of the burn in period
  runtime<-500                  # How many time steps after burn-in?
  times<-seq(0,maxtime1,by=1)   # Vector from 1 to the number of time steps
  mid<-27                       # mean temperature across all reefs at start of simulation.
  range<-3                      # Half of the range of temperatures across reefs at start of simulation (will be added and subtracted from mean to set hottest and coldest temps)
  MPAamount<-0.2                # Proportion of reef network to be managed for local stressors
  cgrow<-1                      # Modifier for coral growth rates - for exploration
  agrow<-1                      # modifier for algae growth rates - for exploration
  tempchange<-3                 # Change in temeprature between historic and future climate conditions
  temp.stoch<-1.2               # Scale of temperature stochasticity
  mdim<-size                    # Number of rows and columns in mvnorm matrix
  lindec<-exp(seq(0,-5,length=mdim))  # Set covariance relationships among sites
  ma<-matrix(0,nrow=mdim,ncol=mdim)   #Create empty matrix to store covariance
  ma[,1]<-lindec                      # Set first column as full range of covariances
  for(i in 2:mdim){                   # Set the lower triangle of all columns in mvnorm matrix                             
    ma[-c(1:(i-1)),i]<-lindec[1:(mdim-(i-1))]
  }
  ma<-ma+t(ma)-diag(nrow(ma))  # Fill in upper triangle
  is.positive.definite(ma)     # Check if positive definite
  sds<-rep(.2*temp.stoch,size) # Scale standard deviations of temperature anomalies
  b<-sds%*%t(sds)              # Variance
  spatialtemp<-b*ma            # Scale the variance covariance matrix for temperature anomalies  
  
  anoms.burn<-matrix(rnorm(burnin,0,temp.stoch),nrow=size,ncol=burnin,byrow=T)+matrix(rmvnorm(burnin,rep(0,size),spatialtemp),nrow=size,byrow=T) # Annual temperature anomalies duringn the burn in period
  anoms.runs<-matrix(rnorm(runtime,0,temp.stoch),nrow=size,ncol=runtime,byrow=T)+matrix(rmvnorm(runtime,rep(0,size),spatialtemp),nrow=size,byrow=T) # annual temperature anomalies during the climate change/management period
  algaemort<-matrix(runif((runtime+burnin)*size,0.05,.3),nrow=size,ncol=runtime+burnin,byrow=T) # Annual algae morality values by reef location across entire time series (these will be over-written if the reef is in a managed area)
  
  #==============================================================================
  # Generate starting conditions and merge into single vector or use in ode.1D()
  #   - Temperatures (specify scenario)
  #   - Densities by species
  #   - Trait (optimum temperatures, specify scenario)
  #
  # Create vector of names for each "State"
  #   - species names, then trait names, then temperature
  #==============================================================================
  temps<-generate.temps(size=size,mid=mid,range=range,temp.scenario="linear") # Generate intial thermal landscape (linear arrangement)
  sppstate<-generate.state(size=size,nsp=nsp,dens=0.25,random=T) # Generate initial species abundances
  traitstate<-generate.traits(nsp=nsp,size=size,mid=mid,temps=temps,range=range,trait.scenario="perfect.adapt") # Generate initial trait distribution (perfectly adapted)
  allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps) # Combine all state variables into a single vector
  
  allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps") # names of state variables
  
  #==============================================================================
  # What types of species (will they be protected by an MPA?)
  #   - Set species that are protected by MPA to sptype=1
  #   - Set species that are NOT protected by MPA to sptype=2
  #   - set MPA value to zero for all sites in burn in, then to 1 for managed sites during experiments
  #==============================================================================
  species<-c("C1","C2","MA")
  sptype<-c(1,1,2)
  mpa<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy="none")
  
  #==============================================================================
  # Define model parameters
  #   - nsp: Number of Species
  #   - mpa: vector of logical values for whether each reef is managed or not
  #   - species: vector of species names
  #   - sptype: vector of species types (1 for coral, 2 for algae)
  #   - V: Genetic variance for each species
  #   - D: Dispersal for each species
  #   - rmax: Maximum growth rate for each species
  #   - alphas: species interaction (competition) matrix (nsp x nsp matrix)
  #   - m: mortality rate for each species
  #   - w: temperature tolerance or plasticity (nsp)
  #   - annual.temp.change: rate of temperature change 
  #   - maxtemp: maximum temperature (sigmoid temp change scenario)
  #   - Nmin: minimum density (for computational happiness)
  #   - deltax: step size in space for derivitives
  #   - timemod: used to offset indices after burn-in
  #==============================================================================
  
  parms<-list(
    nsp=nsp,
    mpa=mpa,
    species=species,
    sptype=sptype,
    V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
    D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
    rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
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
  #     anoms: temperature anomalies
  #     algaemort: algae mortality rates
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

        print("B")  
        
        dspp<-spps   # create storage for change in density state values
        dtraits<-traits # create storage for change in trait state values
        
        for(i in 1:nsp)
        { 
          print("C")
          
          if(i==3) amort<-algaemort[,k+timemod]
          dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                         mort=amort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal") # Calculate change in species cover
          dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                            mort=amort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                            mortality.model="tempvary",spp=species[i],growth="normal") # Calculate change in trait
          
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
  
  out<-coral_trait_stoch(t=maxtime1,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                         temp.change="const",stoch.temp=T,anoms=anoms.burn,algaemort=algaemort)
  
  temps<-out$ts[(2*nsp*size+1):(2*nsp*size+size),maxtime1] # Set final temp state as starting temp state for experiments
  sppstate<-matrix(out$ts[1:(nsp*size),maxtime1],nrow=nsp,ncol=size,byrow=T) # Set final species state from burn in as starting state for experiments
  traitstate<-out$ts[(nsp*size+1):(2*nsp*size),maxtime1]  # set final trait state from burn in as starting trait state for experiments
  allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps) # Concatenate state variables into single vector
  
  allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps") # names for state vector

  strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","random") # Name different management strategies

  nstrat<-length(strategies)   # number of strategies
  mpa<-matrix(NA,nrow=size,ncol=length(strategies)) # Storage matrix for mpa arrangement under different strategies
  
  out2<-array(NA,dim=c(420,500,length(strategies))) # storage array for model output
  
  for(i in 1:length(strategies))
  {
    
    mpa[,i]<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy=strategies[i],priordata=out$ts) # Select sites for management/protection
    maxtime2<-500
    
    #==============================================================================
    # Define model parameters
    #   - nsp: Number of Species
    #   - mpa: vector of logical values for whether each reef is managed or not
    #   - species: vector of species names
    #   - sptype: vector of species types (1 for coral, 2 for algae)
    #   - V: Genetic variance for each species
    #   - D: Dispersal for each species
    #   - rmax: Maximum growth rate for each species
    #   - alphas: species interaction (competition) matrix (nsp x nsp matrix)
    #   - m: mortality rate for each species
    #   - w: temperature tolerance or plasticity (nsp)
    #   - annual.temp.change: rate of temperature change 
    #   - maxtemp: maximum temperature (sigmoid temp change scenario)
    #   - Nmin: minimum density (for computational happiness)
    #   - deltax: step size in space for derivitives
    #   - timemod: used to offset indices after burn-in
    #==============================================================================
    parms<-list(
      nsp=nsp, 
      mpa=mpa[,i], 
      species=species,
      sptype=sptype,
      V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
      D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
      rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),

      alphas=matrix(c(1,1.1,1,
                      1,1,1,
                      .8,.8,1),nrow=nsp,ncol=nsp,byrow=T),
      m=as.numeric(c(t(start.parms[rownames(start.parms)=="M.normal",((1:nsp))]))),
      
      w=c(t(as.numeric(start.parms[rownames(start.parms)=="w",((1:nsp))]))),
      annual.temp.change=.011,
      maxtemp=30,
      Nmin=10^-6,
      deltax=1,
      timemod=burnin
    )
    
    modout<-coral_trait_stoch(t=maxtime2,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "modout"
                              temp.change="sigmoid",stoch.temp=T,anoms=anoms.runs,algaemort=algaemort)
    out2[,1:maxtime2,i]<-modout$ts # save the time series of species cover and reef specific temperatures
    
  }
  
  outputlist<-list(out2,mpa) # Store model output in list form
  
}

registerDoSEQ() # Stops parallel processing
Sys.time()-start.time # Print duration of model run
