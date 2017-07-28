library("RStudioAMI")

passwd()

unlinkDropbox()
linkDropbox()

excludeSyncDropbox("*")
includeSyncDropbox("AWSCoralRuns")

install.packages(c("foreach","doParallel","deSolve","numDeriv","rbenchmark",
                   "doSNOW","fitdistrplus","logitnorm","matrixcalc","mvtnorm",
                   "mail"))


# Loop stochastic Norberg Model across parameter and management scenarios
#=========================================================================
# Read in data files and source functions
source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1) # Matrix of parameter values to be filled from "scenarios"
scenarios<-read.csv("Norberg_StartingParameters_iterationlist_reservesize_SA.csv",header=T)[1:45,] # Each row contains parameter values for one scenario. Column headers identify parameter.
nscen<-nrow(scenarios)
uniqueburns<-length(unique(scenarios$Scenario))
modname<-"Base"

# Set initialization parameters
nsp<-3                        # how many species in model?
size<-60                      # how many reefs in model?
burnin<-1500                  # How long to run burn-in period?
runtime<-500                 # How long to run management/environmental change period?
mid<-27                       # mean temperature across all reefs at start of simulation.
range<-3                      # range of temperatures across reefs at start of simulation
sptype<-c(1,1,2) 
species<-c("C1","C2","MA")
iterations<-50 # Number of stochastic iterations
#monitor.yrs<-c(10,20,30,40,50,100,250,500)
monitor.yrs<-c(20,50,100,500)
coralgrowth<-1
algaegrowth<-1
# MPAamount<-0.2
temp.stoch<-1.2

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

sds<-rep((.2*temp.stoch),size) # SD of intersite temperature variation
b<-sds%*%t(sds)
spatialtemp<-b*ma # Scales MVN matrix to 0.2 * temp.stoch when temp.stoch=1.2
anoms.burn<-array(NA,dim=c(size,burnin,iterations))
anoms.runs<-array(NA,dim=c(size,runtime,iterations))
algaemort<-array(NA,dim=c(size,(runtime+burnin),iterations))
for(i in 1:iterations)
{ 
  set.seed(i)
  
  anoms.burn[,,i]<-matrix(rnorm(burnin,0,temp.stoch),nrow=size,ncol=burnin,byrow=T)+matrix(rmvnorm(burnin,rep(0,size),spatialtemp),nrow=size,byrow=T)
  
  anoms.runs[,,i]<-matrix(rnorm(runtime,0,temp.stoch),nrow=size,ncol=runtime,byrow=T)+matrix(rmvnorm(runtime,rep(0,size),spatialtemp),nrow=size,byrow=T)
  
  algaemort[,,i]<-matrix(runif((runtime+burnin)*size,0.05,.3),nrow=size,ncol=runtime+burnin)
  
  print(i)
}





# Run model iterations in parallel across all cores on computer
cl<-makeCluster(detectCores(),outfile=paste(modname,"_BurnModelOutput_",format(Sys.Date(),format="%m%d%y"),".txt",sep=""))
registerDoParallel(cl)
timestart<-Sys.time()
#tsout<-foreach(qqq=1:nrow(scenarios),.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{
tsout2<-foreach(qqq=1:uniqueburns,.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{

  start.parms<-matrix(scenarios[qqq,][-c(1,(ncol(scenarios)-1),ncol(scenarios))],nrow=5,ncol=3,byrow=F) # fill start.parms matrix from scenarios file
  colnames(start.parms)<-colnames(start.parm.mat) # identify columns
  rownames(start.parms)<-rownames(start.parm.mat) # identify rows
  
  a<-runMod(q=qqq,niter=iterations,burnin=burnin,runtime=runtime,spatialtemp=spatialtemp,
            monitor.yrs=monitor.yrs,burnonly=T,anoms.burn=anoms.burn,anoms.runs=anoms.runs,
            nsp=nsp,size=size,sptype=sptype,species=species,mid=mid,range=range,
            cgrow=coralgrowth,agrow=algaegrowth,amort=algaemort,tscen="linear") # Run the model for given scenario
  a
}

registerDoSEQ() # Close out cluster/parallel operations
stopCluster(cl)
Sys.time()-timestart
setwd("~/Dropbox/AWSCoralRuns")
save.image(paste(modname,"_CoralSimulationsBurnIn",format(Sys.Date(),format="%m%d%y"),".RData",sep=""),compress=T)

# Burn-in: Test Run 9 scenarios x 50 iterations = 450 simulations, takes ~2.5 hrs
strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","portfolioGreedy","random") # Name different management strategies
nstrat<-length(strategies)

#==================
mpaset<-array(0,dim=c(size,nstrat,iterations,nrow(scenarios)))
# 
# for(i in 1:nrow(scenarios))
# {
#   for(j in 1:iterations)
#   {
#     subdat<-tsout2[[scenarios$Scenario[i]]][[j]]$ts
#     for(k in 1:nstrat)
#     {
#       print(scenarios$reserve[i])
#       mpaset[,k,j,i]<-setMPA(temps=subdat[(6*size+1):(7*size),ncol(subdat)],Nall=matrix(subdat[1:(nsp*size),ncol(subdat)],nrow=nsp,ncol=size,byrow=T),
#                              sptype=sptype,size=size,amount=scenarios$reserve[i],strategy=strategies[k],priordata=subdat)
#       
#     # print(paste("i = ", i,"; j = ",j,"; k = ",k))
#     }
#   }
#   
# }


cl<-makeCluster(detectCores(),outfile=paste(modname,"_SetMPAOutput_",format(Sys.Date(),format="%m%d%y"),".txt",sep=""))
registerDoParallel(cl)
timestart<-Sys.time()
mpaset2<-foreach(i=1:nrow(scenarios)) %dopar%{
  mpaset<-mpaset<-array(0,dim=c(size,nstrat,iterations))
  print(i)
  for(j in 1:iterations)
  {
    
    subdat<-tsout2[[scenarios$Scenario[i]]][[j]]$ts
    for(k in 1:nstrat)
    {
      print(scenarios$reserve[i])
      mpaset[,k,j]<-setMPA(temps=subdat[(6*size+1):(7*size),ncol(subdat)],Nall=matrix(subdat[1:(nsp*size),ncol(subdat)],nrow=nsp,ncol=size,byrow=T),
                           sptype=sptype,size=size,amount=scenarios$reserve[i],strategy=strategies[k],priordata=subdat)
      
      # print(paste("i = ", i,"; j = ",j,"; k = ",k))
    }
  }
  
  mpaset
  
}
registerDoSEQ() # Close out cluster/parallel operations
stopCluster(cl)
Sys.time()-timestart

mpaset<-array(0,dim=c(size,nstrat,iterations,nrow(scenarios)))
for(i in 1:nrow(scenarios))
{
  mpaset[,,,i]<-mpaset2[[i]]
}


save(mpaset,file=paste(modname,"_CoralSimulationsMPASetOutData",format(Sys.Date(),format="%m%d%y"),".RData",sep=""))



# Run model iterations in parallel across all cores on computer
cl<-makeCluster(detectCores(),outfile=paste(modname,"_RunModelOutput_",format(Sys.Date(),format="%m%d%y"),".txt",sep=""))
registerDoParallel(cl)
timestart<-Sys.time()
qlist<-rep((seq(1,9)%%9+c(rep(0,8),9)),length(unique(scenarios[,18])))
#tsout<-foreach(qqq=1:nrow(scenarios),.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{
tsoutRuns<-foreach(qqq=1:nrow(scenarios),.packages=c("deSolve","numDeriv","mvtnorm")) %dopar%{
  newdat<-array(NA,dim=c(nrow(tsout2[[1]][[1]]$ts),ncol(tsout2[[1]][[1]]$ts),iterations))
  
  for(i in 1:iterations)
  {
    newdat[,,i]<-tsout2[[qlist[qqq]]][[i]]$ts
  }
  #print("Atest")
  start.parms<-matrix(scenarios[qqq,][-c(1,(ncol(scenarios)-1),ncol(scenarios))],nrow=5,ncol=3,byrow=F) # fill start.parms matrix from scenarios file
  colnames(start.parms)<-colnames(start.parm.mat) # identify columns
  rownames(start.parms)<-rownames(start.parm.mat) # identify rows
  #print("BTEST")
  a<-runMod(q=qqq,niter=iterations,burnin=burnin,runtime=runtime,spatialtemp=spatialtemp,
            monitor.yrs=monitor.yrs,burnonly=F,anoms.burn=anoms.burn,
            anoms.runs=anoms.runs,mpasetuse=mpaset[,,,qqq],priordata=newdat,
            cgrow=coralgrowth,agrow=algaegrowth,nsp=nsp,size=size,sptype=sptype,
            species=species,mid=mid,range=range,amort=algaemort[,-c(1:burnin),],
            tscen="linear") # Run the model for given scenario
  a
}

registerDoSEQ() # Close out cluster/parallel operations
stopCluster(cl)
Sys.time()-timestart

setwd("~/Dropbox/AWSCoralRuns")

#rm(tsout2)
#save.image(paste(modname,"_First16_Scenarios_Runs",format(Sys.Date(),format="%m%d%y"),".RData",sep=""))
#tsoutRunsRandomtemps<-tsoutRuns
setwd("~/")
save(tsoutRuns,file=paste(modname,"_CoralRunsOut",format(Sys.Date(),format="%m%d%y"),".RData",sep=""),compress=T)


#==========================================================================


runMod<-function(q,niter,burnin,runtime,spatialtemp,monitor.yrs,
                 burnonly,anoms.burn,anoms.runs,mpasetuse=NULL,priordata=NULL,agrow=1,cgrow=1,
                 size,nsp,sptype,species,mid,range,amort,tscen="linear")
{ 

  print(q)
  strategies<-c("hot","cold","hotcold","highcoral","lowcoral","space","portfolio","portfolioGreedy","random") # Name different management strategies
  nstrat<-length(strategies)
  
  #stoch.time.series<-matrix(data=NA,nrow=(2*nsp*size+size),ncol=maxtime) # Create storage matrix for time series of state values at each location
#   output.list<-list("C1"=matrix(NA,nrow=nstrat,ncol=(length(monitor.yrs)*2)),"C2"=matrix(NA,nrow=nstrat,ncol=(length(monitor.yrs)*2)),
#                     "Call"=matrix(NA,nrow=nstrat,ncol=(length(monitor.yrs)*2)),"MA"=matrix(NA,nrow=nstrat,ncol=(length(monitor.yrs)*2)),"Shannon"=matrix(NA,nrow=nstrat,ncol=(length(monitor.yrs)*2)))
    #print("A")
  if(burnonly==T)
  {
    #   set.seed(seeds[q])
    #set.seed(iters)
    #set.seed(1)
    maxtime1<-burnin              # how many time steps?
    
    times<-seq(0,maxtime1,by=1)   # Vector from 1 to the number of time steps
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
    temps<-matrix(generate.temps(size=size,mid=mid,range=range,temp.scenario=tscen),ncol=niter,nrow=size,byrow=F) 
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
    
    #allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")
    #print("C")
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
      spatialtemp=spatialtemp
      
    )
    
    #print("D")
    start.time<-Sys.time()
    out1<-apply(allstate,MARGIN=2,FUN=coral_trait_stoch.vec,t=maxtime1,parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                temp.change="const",stoch.temp=T,burn=T,anoms=anoms.burn,amort=amort)
    Sys.time()-start.time  
    return(out1)
    
  }
    
    
    #print("BurnCheck")
    #print(out$ts[1,maxtime])
  if(burnonly==F)
  {
    tscen<-NA
    print("RunCheck")
    maxtime<-runtime                   # how many time steps?
    times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
    output.array<-array(NA,dim=c(dim(priordata)[1],14*length(monitor.yrs)+3,niter,nstrat))
    
    #temps<-generate.temps(size=size,mid=mid,range=range,temp.scenario="linear") 
    temps<-priordata[(2*nsp*size+1):((2*nsp+1)*size),dim(priordata)[2],]
    #print(temps)
    #sppstate<-out1$ts[(1:(nsp*size)),maxtime1]
    sppstate<-priordata[1:(nsp*size),dim(priordata)[2],]
    print("SPPCheck")
    #traitstate<-out1$ts[(nsp*size+1):(2*nsp*size),maxtime1]
    traitstate<-priordata[(nsp*size+1):(2*nsp*size),dim(priordata)[2],]
    print("TraitCheck")
    allstate1<-rbind(sppstate,traitstate,temps)
    print("RBINDCheck")
    #allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")
    
    for(strategy in 1:nstrat)
    {
      mpa<-mpasetuse[,strategy,]
      print(strategy)
      #print(dim(allstate1))
      print("MPACheck")
      allstate<-rbind(allstate1,mpa,seq(1,niter))
     # print(dim(allstate))
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
        spatialtemp=spatialtemp
      )
      #print(allstate[,1])
      print("E")
#       out<-coral_trait_stoch(t=maxtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
#                              temp.change="sigmoid",stoch.temp=T,burn=F,anoms=anoms.runs)
      print("yo")
      out<-apply(allstate,MARGIN=2,FUN=coral_trait_stoch.vec,t=maxtime,parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                  temp.change="sigmoid",stoch.temp=T,burn=F,anoms=anoms.runs,amort=amort)
      for(iter in 1:niter)
      {
        print("a")
        output.array[,,iter,strategy]<-out[[iter]]$ts
        print("b")
      }
      #if(strategy==1) output.list$ExampleTimeSeries[,]<-out$ts
      
      #       c1means<-colMeans(out$ts[1:size,])[monitor.yrs]
      #       c2means<-colMeans(out$ts[(1+size):(2*size),])[monitor.yrs]
      #       callmeans<-c1means+c2means
      #       mameans<-colMeans(out$ts[(2*size+1):(3*size),])[monitor.yrs]
      #       
      #       c1sd<-apply(out$ts[1:size,],MARGIN=2,FUN=sd)[monitor.yrs]
      #       c2sd<-apply(out$ts[(1+size):(2*size),],MARGIN=2,FUN=sd)[monitor.yrs]
      #       callsd<-apply((out$ts[1:(size),]+out$ts[(size+1):(2*size),]),MARGIN=2,FUN=sd)[monitor.yrs]
      #       masd<-apply(out$ts[(2*size+1):(3*size),],MARGIN=2,FUN=sd)[monitor.yrs]
      #       
      #       output.list$C1[strategy,seq(1,2*length(monitor.yrs)-1,by=2)]<-c1means
      #       output.list$C2[strategy,seq(1,2*length(monitor.yrs)-1,by=2)]<-c2means
      #       output.list$Call[strategy,seq(1,2*length(monitor.yrs)-1,by=2)]<-callmeans
      #       output.list$MA[strategy,seq(1,2*length(monitor.yrs)-1,by=2)]<-mameans
      #       
      #       output.list$C1[strategy,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-c1sd
      #       output.list$C2[strategy,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-c2sd
      #       output.list$Call[strategy,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-callsd
      #       output.list$MA[strategy,(seq(1,2*length(monitor.yrs)-1,by=2)+1)]<-masd
      
    }
    
    
    return(output.array)
    #return(extinct.array)
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
                                stoch.temp=TRUE,burn=TRUE,anoms,amort){
  
  with(as.list(parms), {  # Load parameter values
    if(burn==TRUE){
      anoms<-anoms
      time.series<-matrix(NA,nrow=(2*nsp*size+size),ncol=t)
      time.series[,1]<-y[1:(2*nsp*size+size)]  # set initial state
      mpa<-y[(7*size+1):(8*size)]
      stoch.it<-y[length(y)]
      for(k in 2:t)
      {
        # print(k) # Print counter
        
        
        spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
        traits<-matrix(time.series[(size*nsp+1):(length(time.series[,1])-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
        temps<-time.series[(length(time.series[,1])-(size-1)):length(time.series[,1]),k-1]+stoch.temp*anoms[,k-1,stoch.it] # Extract temperature state values
        #print("Anomalies Check")
        #       print(temps)
        #       print(anoms)
        dspp<-spps   # create storage for change in density state values
        dtraits<-traits # create storage for change in trait state values
        
        
        for(i in 1:nsp)
        { 
          #print("AA")
          if(i==3) almort<-amort[,k,stoch.it]
          else almort<-rep(0,size)
          dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                         mort=almort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
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
        time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
        
      }
      return(list("ts"=time.series[,(t/2+1):t]))
      
    }
    
    if(burn==FALSE)
    {
      time.series.monitor<-matrix(NA,nrow=(2*nsp*size+size),ncol=(14*length(monitor.yrs)+3))# Create storage matrix for time series of state values at each location
      #anoms<-matrix(NA,nrow=size,ncol=(t-1))
      
      time.series<-y[1:(2*nsp*size+size)]  # set initial state
      time.series.monitor[,1]<-time.series
      mpa<-y[(7*size+1):(8*size)]
     # print(mpa)
      stoch.it<-y[length(y)]
      #print(stoch.it)
      extinctions.001<-matrix(0,nrow=3,ncol=(t-1))
      extinctions.01<-matrix(0,nrow=3,ncol=(t-1))
      extinctions.05<-matrix(0,nrow=3,ncol=(t-1))
      extinctions2.001<-matrix(0,nrow=3,ncol=(t-1))
      extinctions2.01<-matrix(0,nrow=3,ncol=(t-1))
      extinctions2.05<-matrix(0,nrow=3,ncol=(t-1))
      for(k in 2:t)
      {
       # print(k) # Print counter
        
        
        spps<-matrix(time.series[(1:(size*nsp))],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
        #print(dim(spps))
        traits<-matrix(time.series[(size*nsp+1):(length(time.series)-size)],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
       # print(dim(traits))
        temps<-time.series[(length(time.series)-(size-1)):length(time.series)]
        #print(length(temps))
        current.anom<-stoch.temp*anoms[,k-1,stoch.it]
        #print(length(current.anom))
        temps<-temps+current.anom # Extract temperature state values
        #print("Anomalies Check")
        #       print(temps)
        #       print(anoms)
        dspp<-spps   # create storage for change in density state values
        dtraits<-traits # create storage for change in trait state values
        #print(spps)
        #print(traits)
        for(i in 1:nsp)
        { 
          #print("AA")
          if(i==3) almort<-amort[,k,stoch.it]
          else almort<-rep(0,size)
          dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],
                         mort=almort,TC=temps,alphas=alphas[i,],delx=1,mpa=mpa,mortality.model="tempvary",spp=species[i],growth="normal")
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
        
        time.series<-time.series+dsp # Store new state value
        extinctions.001[,k-1]<-extinctrisk(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),size=size,threshold=0.001)
        extinctions.01[,k-1]<-extinctrisk(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),size=size,threshold = 0.01)
        extinctions.05[,k-1]<-extinctrisk(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),size=size,threshold=0.05)
        
        extinctions2.001[,k-1]<-extinctrisk2(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),threshold=0.001)
        extinctions2.01[,k-1]<-extinctrisk2(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),threshold = 0.01)
        extinctions2.05[,k-1]<-extinctrisk2(Nall=matrix(time.series[1:(3*size)],nrow=3,ncol=size,byrow=T),threshold=0.05)
        
        print(extinctions.05[,k-1])
        if(k %in% monitor.yrs) {
          time.series.monitor[,which(monitor.yrs==k)+1]<-time.series
          time.series.monitor[,c((length(monitor.yrs)*1)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.001[,1:(k-1)])/(k-1))[1]
          time.series.monitor[,c((length(monitor.yrs)*2)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.001[,1:(k-1)])/(k-1))[2]
          print("Check1")
          time.series.monitor[,c((length(monitor.yrs)*3)+(which(monitor.yrs==k))+1)]<-(sum(colSums(extinctions.001[1:2,1:(k-1)])==2)/(k-1))
          
          time.series.monitor[,c((length(monitor.yrs)*4)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.01[,1:(k-1)])/(k-1))[1]
          time.series.monitor[,c((length(monitor.yrs)*5)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.01[,1:(k-1)])/(k-1))[2]
          print("Check1")
          time.series.monitor[,c((length(monitor.yrs)*6)+(which(monitor.yrs==k))+1)]<-(sum(colSums(extinctions.01[1:2,1:(k-1)])==2)/(k-1))
          
          time.series.monitor[,c((length(monitor.yrs)*7)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.05[,1:(k-1)])/(k-1))[1]
          time.series.monitor[,c((length(monitor.yrs)*8)+(which(monitor.yrs==k))+1)]<-(rowSums(extinctions.05[,1:(k-1)])/(k-1))[2]
          print("Check1")
          time.series.monitor[,c((length(monitor.yrs)*9)+(which(monitor.yrs==k))+1)]<-(sum(colSums(extinctions.05[1:2,1:(k-1)])==2)/(k-1))
          print("OK")
        }

        #print(extinctions[,k-1])
      }
      #print("AAAAA")
      #print(extinctions[1,])
      #print(1 %in% extinctions[1,])
      #print(which(extinctions[1,]==1,arr.ind=T)[1])
      time.series.monitor[,10*length(monitor.yrs)+2]<-ifelse(sum(extinctions.001[1,])>0,which(extinctions.001[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+3]<-ifelse(sum(extinctions.001[2,])>0,which(extinctions.001[2,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+5]<-ifelse(sum(extinctions.01[1,])>0,which(extinctions.01[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+6]<-ifelse(sum(extinctions.01[2,])>0,which(extinctions.01[2,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+8]<-ifelse(sum(extinctions.05[1,])>0,which(extinctions.05[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+9]<-ifelse(sum(extinctions.05[2,])>0,which(extinctions.05[2,]==1)[1],NA)
      #print("check2")
      both.extinct.001<-which(colSums(extinctions.001[1:2,])==2)
      both.extinct.01<-which(colSums(extinctions.01[1:2,])==2)
      both.extinct.05<-which(colSums(extinctions.05[1:2,])==2)
      # time.series.monitor[,10*length(monitor.yrs)+4]<-ifelse(sum(both.extinct.001[1,])>0,which(both.extinct.001[1,]==1)[1],NA)
      # time.series.monitor[,10*length(monitor.yrs)+7]<-ifelse(sum(both.extinct.01[1,])>0,which(both.extinct.01[1,]==1)[1],NA)
      # time.series.monitor[,10*length(monitor.yrs)+10]<-ifelse(sum(both.extinct.05[1,])>0,which(both.extinct.05[1,]==1)[1],NA)
      
      time.series.monitor[,10*length(monitor.yrs)+11]<-ifelse(sum(extinctions2.001[1,])>0,which(extinctions2.001[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+12]<-ifelse(sum(extinctions2.001[2,])>0,which(extinctions2.001[2,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+14]<-ifelse(sum(extinctions2.01[1,])>0,which(extinctions2.01[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+15]<-ifelse(sum(extinctions2.01[2,])>0,which(extinctions2.01[2,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+17]<-ifelse(sum(extinctions2.05[1,])>0,which(extinctions2.05[1,]==1)[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+18]<-ifelse(sum(extinctions2.05[2,])>0,which(extinctions2.05[2,]==1)[1],NA)
      #print("check2")
      both.extinct2.001<-which(colSums(extinctions2.001[1:2,])==2)
      both.extinct2.01<-which(colSums(extinctions2.01[1:2,])==2)
      both.extinct2.05<-which(colSums(extinctions2.05[1:2,])==2)
      
      #print("OK")
      #print(both.extinct)
      #print(both.extinct[1])
      #print("check3")
      #print(!is.na(both.extinct[1]))
      time.series.monitor[,10*length(monitor.yrs)+13]<-ifelse(length(both.extinct2.001[1])!=0,both.extinct2.001[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+16]<-ifelse(length(both.extinct2.01[1])!=0,both.extinct2.01[1],NA)
      time.series.monitor[,10*length(monitor.yrs)+19]<-ifelse(length(both.extinct2.05[1])!=0,both.extinct2.05[1],NA)
      #print("OK")
      return(list("ts"=time.series.monitor))
      
    }
    
  })
  
  
  
}



