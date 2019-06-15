install.packages(c("foreach","doParallel","deSolve","numDeriv","rbenchmark",
                   "doSNOW","fitdistrplus","logitnorm","matrixcalc","mvtnorm",
                   "mail"))


# Loop stochastic Norberg Model across parameter and management scenarios
#=========================================================================
# Read in data files and source functions
source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1) # Matrix of parameter values to be filled from "scenarios"
scenarios<-read.csv("Norberg_StartingParameters_iterationlist.csv",header=T) # Each row contains parameter values for one scenario. Column headers identify parameter.


# Set initialization parameters
nsp<-3                        # how many species in model?
size<-100                      # how many reefs in model?
burnin<-1500                  # How long to run burn-in period?
runtime<-500                 # How long to run management/environmental change period?
mid<-27                       # mean temperature across all reefs at start of simulation.
range<-5                      # range of temperatures across reefs at start of simulation
sptype<-c(1,1,2) 
species<-c("C1","C2","MA")
iterations<-100 # Number of stochastic iterations
monitor.yrs<-c(20,50,100,500)
coralgrowth<-1
algaegrowth<-1
MPAamount<-0.2

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
anoms.burn<-array(NA,dim=c(size,burnin,50))
anoms.runs<-array(NA,dim=c(size,runtime,50))

for(i in 1:50)
{ 
  set.seed(i)
  
  anoms.burn[,,i]<-matrix(rnorm(burnin,0,1.2),nrow=size,ncol=burnin,byrow=T)+matrix(rmvnorm(burnin,rep(0,size),spatialtemp),nrow=size,byrow=T)
  
  anoms.runs[,,i]<-matrix(rnorm(runtime,0,1.2),nrow=size,ncol=runtime,byrow=T)+matrix(rmvnorm(runtime,rep(0,size),spatialtemp),nrow=size,byrow=T)
  
  print(i)
}


# Run model iterations in parallel across all cores on computer
cl<-makeCluster(9,outfile="BurnModelOutput110116.txt")
registerDoParallel(cl)
timestart<-Sys.time()
tsout2<-foreach(qqq=1:9,.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{

  start.parms<-matrix(scenarios[qqq,][-c(1,ncol(scenarios))],nrow=5,ncol=3,byrow=F) # fill start.parms matrix from scenarios file
  colnames(start.parms)<-colnames(start.parm.mat) # identify columns
  rownames(start.parms)<-rownames(start.parm.mat) # identify rows
  
  a<-runMod(q=qqq,niter=iterations,burnin=burnin,runtime=runtime,spatialtemp=spatialtemp,
            monitor.yrs=monitor.yrs,burnonly=T,anoms.burn=anoms.burn,anoms.runs=anoms.runs,
            nsp=nsp,size=size,sptype=sptype,species=species,mid=mid,range=range,
            cgrow=coralgrowth,agrow=algaegrowth) # Run the model for given scenario
  a
}

registerDoSEQ() # Close out cluster/parallel operations
stopCluster(cl)
Sys.time()-timestart
save.image(paste("CoralSimulationsBurnIn",format(Sys.Date(),format="%m%d%y"),".RData",sep=""))

# Burn-in: Test Run 9 scenarios x 50 iterations = 450 simulations, takes ~2.5 hrs
strategies<-c("hot","cold","highcoral","lowcoral","portfolio","random","none") # Name different management strategies
nstrat<-length(strategies)

#==================
mpaset<-array(0,dim=c(size,7,iterations,9))

for(i in 1:9)
{
  for(j in 1:iterations)
  {
    subdat<-tsout2[[i]][[j]]$ts
    for(k in 1:nstrat)
    {
      mpaset[,k,j,i]<-setMPA(temps=subdat[(6*size+1):(7*size),ncol(subdat)],Nall=matrix(subdat[1:(2*size),ncol(subdat)],nrow=nsp,ncol=size,byrow=T),
                             sptype=sptype,size=size,amount=MPAamount,strategy=strategies[k],priordata=subdat)
      
      print(paste("i = ", i,"; j = ",j,"; k = ",k))
    }
  }
  
}


# Run model iterations in parallel across all cores on computer
cl<-makeCluster(9,outfile="RunModelOutput110116.txt")
registerDoParallel(cl)
timestart<-Sys.time()
tsoutRuns<-foreach(qqq=1:9,.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{
  newdat<-array(NA,dim=c(nrow(tsout2[[1]][[1]]$ts),ncol(tsout2[[1]][[1]]$ts),iterations))
  
  for(i in 1:iterations)
  {
    newdat[,,i]<-tsout2[[qqq]][[i]]$ts
  }
  #print("Atest")
  start.parms<-matrix(scenarios[qqq,][-c(1,ncol(scenarios))],nrow=5,ncol=3,byrow=F) # fill start.parms matrix from scenarios file
  colnames(start.parms)<-colnames(start.parm.mat) # identify columns
  rownames(start.parms)<-rownames(start.parm.mat) # identify rows
  #print("BTEST")
  a<-runMod(q=qqq,niter=iterations,burnin=burnin,runtime=runtime,spatialtemp=spatialtemp,
            monitor.yrs=monitor.yrs,burnonly=F,anoms.burn=anoms.burn,
            anoms.runs=anoms.runs,mpasetuse=mpaset[,,,qqq],priordata=newdat,
            cgrow=coralgrowth,agrow=algaegrowth,nsp=nsp,size=size,sotype=sptype,species=species,mid=mid,range=range) # Run the model for given scenario
  a
}

registerDoSEQ() # Close out cluster/parallel operations
stopCluster(cl)
Sys.time()-timestart
rm(tsout2)
save.image(paste("CoralSimulationsRuns",format(Sys.Date(),format="%m%d%y"),".RData",sep=""))




#==========================================================================


runMod<-function(q,niter,burnin,runtime,spatialtemp,monitor.yrs,
                 burnonly,anoms.burn,anoms.runs,mpasetuse=NULL,priordata=NULL,agrow=1,cgrow=1,
                 size,nsp,sptype,species,mid,range)
{ 


  strategies<-c("hot","cold","highcoral","lowcoral","portfolio","random","none") # Name different management strategies
  nstrat<-length(strategies)      

  if(burnonly==T)
  {

    maxtime1<-burnin              # how many time steps?
    times<-seq(0,maxtime1,by=1)   # Vector from 1 to the number of time steps

    #==============================================================================
    # Generate starting conditions and merge into single vector or use in ode.1D()
    #   - Temperatures (specify scenario)
    #   - Densities by species
    #   - Trait (optimum temperatures, specify scenario)
    #
    # Create vector of names for each "State"
    #   - species names, then trait names, then temperature
    #==============================================================================
    temps<-matrix(generate.temps(size=size,mid=mid,range=range,temp.scenario="linear"),ncol=niter,nrow=size,byrow=F) 
    sppstate<-matrix(NA,nrow=nsp*size,ncol=niter)
    for(i in 1:niter)
    {
      set.seed(i)
      sppstate[,i]<-generate.state(size=size,nsp=nsp,dens=0.25,random=T)
    }
    
    traitstate<-matrix(NA,nrow=nsp*size,ncol=niter)
    for(i in 1:niter)
    {
      traitstate[,i]<-generate.traits(nsp=nsp,size=size,mid=mid,temps=temps[,i],range=range,trait.scenario="perfect.adapt")
    }
    allstate<-rbind(sppstate,traitstate,temps)
    

    #==============================================================================
    # What types of species (will they be protected by an MPA?)
    #   - Set species that are protected by MPA to sptype=1
    #   - Set species that are NOT protected by MPA to sptype=2
    #   - For species where sptype=2, set mpa to a value >1
    #   - For species where sptype=1, set mpa<=1
    #   - mpa will be multiplied by mortality rate
    #==============================================================================

    
    mpa<-setMPA(temps=temps,Nall=matrix(sppstate,nrow=nsp,ncol=size,byrow=T),sptype=sptype,size=size,amount=0.2,
                strategy="none")
    allstate<-rbind(allstate,matrix(mpa,nrow=size,ncol=niter,byrow=F),seq(1,niter))
    
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

      monitor.yrs<-monitor.yrs,
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
      temp.stoch=1, 
      annual.temp.change=.011,
      maxtemp=30,
      Nmin=10^-6,
      deltax=1,
      spatialtemp=spatialtemp
      
    )
    
    start.time<-Sys.time()
    out1<-apply(allstate,MARGIN=2,FUN=coral_trait_stoch.vec,t=maxtime1,parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                temp.change="const",stoch.temp=T,burn=T,anoms=anoms.burn)
    Sys.time()-start.time  
    return(out1)
    
  }
    

  if(burnonly==F)
  {

    maxtime<-runtime                   # how many time steps?
    times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
    output.array<-array(NA,dim=c(dim(priordata)[1],length(monitor.yrs),niter,nstrat))
     
    temps<-priordata[(2*nsp*size+1):((2*nsp+1)*size),dim(priordata)[2],]

    sppstate<-priordata[1:(nsp*size),dim(priordata)[2],]
   
    traitstate<-priordata[(nsp*size+1):(2*nsp*size),dim(priordata)[2],]

    allstate1<-rbind(sppstate,traitstate,temps)

    for(strategy in 1:nstrat)
    {
      mpa<-mpasetuse[,strategy,]

      allstate<-rbind(allstate1,mpa,seq(1,niter))

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

        monitor.yrs<-monitor.yrs,
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
        temp.stoch=1, 
        annual.temp.change=.011,
        maxtemp=30,
        Nmin=10^-6,
        deltax=1,
        spatialtemp=spatialtemp
      )
  
      out<-apply(allstate,MARGIN=2,FUN=coral_trait_stoch.vec,t=maxtime,parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                  temp.change="sigmoid",stoch.temp=T,burn=F,anoms=anoms.runs)

      for(iter in 1:niter)
      {
        output.array[,,iter,strategy]<-out[[iter]]$ts
      }
      
      
    }
    
    
    return(output.array)
    
  }
  
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

coral_trait_stoch.vec<-function(y,t, parms,size,nsp,temp.change=c("const","linear","sigmoid"),
                                stoch.temp=TRUE,burn=TRUE,anoms){
  
  with(as.list(parms), {  # Load parameter values
    if(burn==TRUE){
      anoms<-anoms
      time.series<-matrix(NA,nrow=(2*nsp*size+size),ncol=t)
      time.series[,1]<-y[2*nsp*size+size]  # set initial state
      mpa<-y[(7*size+1):(8*size)]
      stoch.it<-y[length(y)]
      for(k in 2:t)
      {
       
        
        
        spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
        traits<-matrix(time.series[(size*nsp+1):(length(time.series[,1])-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
        temps<-time.series[(length(time.series[,1])-(size-1)):length(time.series[,1]),k-1]+stoch.temp*anoms[,k-1,stoch.it] # Extract temperature state values
       
        dspp<-spps   # create storage for change in density state values
        dtraits<-traits # create storage for change in trait state values
        
        
        for(i in 1:nsp)
        { 
          
          dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                         mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
          
          dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                            mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                            mortality.model="tempvary",spp=species[i],growth="normal")
          
        }
        
        # Calculate change in temperatures across reef
        dtemps<-switch(temp.change,
                       const=rep(0,size),
                       linear=rep(annual.temp.change,size),
                       sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))
        
        dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector        
        time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
        
      }
      return(list("ts"=time.series[,(t/2+1):t]))
      
    }
    
    if(burn==FALSE)
    {
      time.series.monitor<-matrix(NA,nrow=(2*nsp*size+size),ncol=length(monitor.yrs))# Create storage matrix for time series of state values at each location
      
      time.series<-y[1:(2*nsp*size+size)]  # set initial state
      mpa<-y[(7*size+1):(8*size)]
     
      stoch.it<-y[length(y)]
      
      for(k in 2:t)
      {
       
        
        
        spps<-matrix(time.series[(1:(size*nsp))],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
       
        traits<-matrix(time.series[(size*nsp+1):(length(time.series)-size)],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
        temps<-time.series[(length(time.series)-(size-1)):length(time.series)]
       
        current.anom<-stoch.temp*anoms[,k-1,stoch.it]
       
        temps<-temps+current.anom # Extract temperature state values
       
        dspp<-spps   # create storage for change in density state values
        dtraits<-traits # create storage for change in trait state values
        
        for(i in 1:nsp)
        { 
          
          dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                         mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
          
          dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                            mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,Nmin=Nmin,
                            mortality.model="tempvary",spp=species[i],growth="normal")
          
        }
        
        # Calculate change in temperatures across reef
        dtemps<-switch(temp.change,
                       const=rep(0,size),
                       linear=rep(annual.temp.change,size),
                       sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))
        
        dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps) # concatenate changes in state back to single vector        
        
        time.series<-time.series+dsp  # Store new state value
        if(k %in% monitor.yrs) time.series.monitor[,which(monitor.yrs==k)]<-time.series
      }
      return(list("ts"=time.series.monitor))
      
    }
    
  })
  
  
  
}



