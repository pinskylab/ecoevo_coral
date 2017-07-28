source("Norberg_Functions.r")
start.parm.mat<-read.csv("Norberg_StartingParameters_test.csv",header=T,row.names=1)
scenarios<-read.csv("Norberg_StartingParameters_iterationlist.csv",header=T)
#start.parms<-start.parm.mat
start.parms<-matrix(scenarios[5,-1],nrow=5,ncol=3,byrow=F)
colnames(start.parms)<-colnames(start.parm.mat)
rownames(start.parms)<-rownames(start.parm.mat)

library(mvtnorm)


# Making A Change
set.seed(5)
nsp<-3                        # how many species in model?
size<-60                      # how many reefs in model?
maxtime<-1500 
burnin<-maxtime
runtime<-500# how many time steps?
times<-seq(0,maxtime,by=1)    #Vector from 1 to the number of time steps
mid<-28                       # mean temperature across all reefs at start of simulation.
range<-5                      # range of temperatures across reefs at start of simulation
MPAamount<-0.2
cgrow<-1
agrow<-1
tempchange<-2

mdim<-size
lindec<-exp(seq(0,-6,length=mdim))#rev(seq(-0.9,1,length=mdim))
ma<-matrix(0,nrow=mdim,ncol=mdim)
ma[,1]<-lindec
for(i in 2:mdim){
  ma[-c(1:(i-1)),i]<-lindec[1:(mdim-(i-1))]
}
ma<-ma+t(ma)-diag(nrow(ma))
is.positive.definite(ma)


sds<-rep(.24,size)
b<-sds%*%t(sds)
spatialtemp<-b*ma

anoms.burn<-matrix(rnorm(burnin,0,1.2),nrow=size,ncol=burnin,byrow=T)+matrix(rmvnorm(burnin,rep(0,size),spatialtemp),nrow=size,byrow=T)

anoms.runs<-matrix(rnorm(runtime,0,1.2),nrow=size,ncol=runtime,byrow=T)+matrix(rmvnorm(runtime,rep(0,size),spatialtemp),nrow=size,byrow=T)
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
  mpa=mpa,
  V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
  D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
  rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
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
  maxtemp=mid+tempchange,
  Nmin=10^-6,
  deltax=1,
  spatialtemp=spatialtemp
  
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
                            stoch.temp=TRUE,anoms){
  
  with(as.list(parms), {  # Load parameter values
    time.series<-matrix(NA,nrow=(2*nsp*size+size),ncol=t) # Create storage matrix for time series of state values at each location
    dspout<-time.series # Create storage matrix for time series of change in state values at each location
    time.series[,1]<-y  # set initial state
    
    for(k in 2:t)
    {
      print(k) # Print counter
      
      spps<-matrix(time.series[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)   # Extract species density state values
      traits<-matrix(time.series[(size*nsp+1):(length(y)-size),k-1],nrow=nsp,ncol=size,byrow=T)  # Extract trait state values
      temps<-time.series[(length(y)-(size-1)):length(y),k-1]+anoms[,k-1] # Extract temperature state values
      #print(temps)  
#       pcat<-rbinom(1,1,pcatastrophe)           # Binomial draw for annual catastrophe occurance
#       m<-(1-pcat)*m.normal+pcat*m.catastrophe  # set annual mortality rate dependent on presence of catastrophe
#       
      m<-m
      dspp<-spps   # create storage for change in density state values
      dtraits<-traits # create storage for change in trait state values
      #print(w)
      for(i in 1:nsp)
      { 
        #print(w[i]^2)
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
      

      dspout[,k]<-dsp  # Store change in state value
      time.series[,k]<-time.series[,k-1]+dsp  # Store new state value
      
    }
    return(list("ts"=time.series,"dsp"=dspout))
  })


  
}
start.time<-Sys.time()
out<-coral_trait_stoch(t=maxtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                       temp.change="const",stoch.temp=T,anoms=anoms.burn)
Sys.time()-start.time

#   }
# }

temps<-out$ts[(2*nsp*size+1):(2*nsp*size+size),maxtime] 
sppstate<-matrix(out$ts[1:(nsp*size),maxtime],nrow=nsp,ncol=size,byrow=T) 
traitstate<-out$ts[(nsp*size+1):(2*nsp*size),maxtime] 
allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)

allnames<-c(paste("spp",seq(1,nsp),sep=""),paste("opt",seq(1,nsp),sep=""),"temps")
cols<-colorRampPalette(c("red","yellow","blue"))
cols2<-c("dodgerblue","darkorange")
strategies<-c("hot","cold","highcoral","lowcoral","portfolio","random","none") # Name different management strategies
nstrat<-length(strategies)    
mpa<-matrix(NA,nrow=size,ncol=7)
out2<-array(NA,dim=c(420,500,7))
for(i in 1:7)
{
  mpa[,i]<-setMPA(temps=temps,Nall=sppstate,sptype=sptype,size=size,amount=MPAamount,strategy=strategies[i],priordata=out$ts)
  maxtime2<-500
  
  parms<-list(
    nsp=nsp,
    mpa=mpa[,i],
    V=as.numeric(c(t(start.parms[rownames(start.parms)=="V",((1:nsp))]))),
    D=as.numeric(c(t(start.parms[rownames(start.parms)=="D",((1:nsp))]))),
    rmax=c(cgrow,cgrow,agrow)*as.numeric(c(t(start.parms[rownames(start.parms)=="Rmax",((1:nsp))]))),
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
    maxtemp=mid+tempchange,
    Nmin=10^-6,
    deltax=1,
    spatialtemp=spatialtemp
    
  )
  out2[,,i]<-coral_trait_stoch(t=runtime,y=allstate, parms=parms,size=size,nsp=nsp,     # Runs the model and stores output as "out"
                         temp.change="sigmoid",stoch.temp=T,anoms=anoms.runs)$ts
  
  j<-1 # start counter
#   while(j<=max(hard.corals)) # loop through the hard coral species and plot the average density across the reef through time
#   {
#                          
#     lines(seq(1,runtime)+1500,(colSums(out$ts[((j-1)*size+1):(j*size),])/size),col=cols2[j],lwd=1)
#     j<-j+1
#   }
# lines(seq(1,runtime)+1500,(colSums(out$ts[((1-1)*size+1):(max(hard.corals)*size),])/size),col=cols(7)[i],lwd=2) # Plot the total hard coral density across the reef through time
  # 
}


hard.corals<-which(sptype==1) # Which of the species are hard corals?
macro.algae<-which(sptype==2) # Which of the species are macroalgae?

windows()  # opens graphics window
layout(matrix(c(1,
                2,
                3,
                4),nrow=4,ncol=1))  # Layout panels according to nsp
par(mar=c(3,3,0,0))
j<-1 # start counter

plot(seq(1,maxtime),(colSums(out$ts[((1-1)*size+1):(max(hard.corals)*size),])/size),ylim=c(0,1),xlim=c(0,2000),
     col="black",type="l",lwd=2) # Plot the total hard coral density across the reef through time
   
for(i in 1:7)
{
  
  lines(seq(1,runtime)+1500,(colSums(out2[((1-1)*size+1):(max(hard.corals)*size),,i])/size),col=cols(7)[i],lwd=2) # Plot the total hard coral density across the reef through time
  
}
legend("top",horiz=T,legend=strategies,lty=1,lwd=2,col=cols(7),bty="n",cex=2)


plot(seq(1,maxtime),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1,2000),type="l",lwd=2)
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:1500)
{
  
  
  
  weight.temps[,j]<-out$ts[(3*size+1):(4*size),j]*(out$ts[(1*size+1):(2*size),j]/(colSums(out$ts[(1*size+1):(2*size),])[j]))
  
  
}

lines(seq(1,maxtime),colSums(weight.temps),col="darkorange",lty=2,lwd=2)

#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2)
for(i in 1:7)
{
  weight.temps<-out2[(1*size+1):(2*size),,i]
  for(j in 1:500)
  {
    
    

      weight.temps[,j]<-out2[(3*size+1):(4*size),j,i]*(out2[(1*size+1):(2*size),j,i]/(colSums(out2[(1*size+1):(2*size),,i])[j]))

    
  }
  lines(seq(1,runtime)+1500,colMeans(out2[(6*size+1):(7*size),,i]),xlim=c(1,2000),lwd=2)
#   lines(seq(1,runtime)+1500,colMeans(out2[(3*size+1):(4*size),,i]),col=cols(7)[i],lty=2,lwd=2)
 # lines(seq(1,runtime)+1500,colMeans(out$ts[(4*size+1):(5*size),,i]),col="dodgerblue",lty=2)
 lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols(7)[i],lty=2,lwd=2)
  
}


plot(seq(1,maxtime),colMeans(out$ts[(6*size+1):(7*size),]),ylim=c(25,35),xlim=c(1,2000),type="l",lwd=2)
#lines(seq(1,maxtime),colMeans(out$ts[(3*size+1):(4*size),]),col="darkorange",lty=2)
#lines(seq(1,maxtime),colMeans(out$ts[(4*size+1):(5*size),]),col="dodgerblue",lty=2,lwd=2)
weight.temps<-out$ts[(1*size+1):(2*size),]
for(j in 1:1500)
{
  
  
  
  weight.temps[,j]<-out$ts[(4*size+1):(5*size),j]*(out$ts[(2*size+1):(3*size),j]/(colSums(out$ts[(2*size+1):(3*size),])[j]))
  
  
}



lines(seq(1,maxtime),colSums(weight.temps),col="dodgerblue",lty=2,lwd=2)

for(i in 1:7)
{
  weight.temps<-out2[(4*size+1):(5*size),,i]
  for(j in 1:500)
  {

      weight.temps[,j]<-out2[(4*size+1):(5*size),j,i]*(out2[(2*size+1):(3*size),j,i]/(colSums(out2[(2*size+1):(3*size),,i])[j]))

  
  }
  lines(seq(1,runtime)+1500,colMeans(out2[(6*size+1):(7*size),,i]),ylim=c(25,30),xlim=c(1,2000),lwd=2)
  #lines(seq(1,runtime)+1500,colMeans(out$ts[(3*size+1):(4*size),,i]),col="darkorange",lty=2)
  #lines(seq(1,runtime)+1500,colMeans(out2[(4*size+1):(5*size),,i]),col=cols(7)[i],lty=2,lwd=2)
  
  lines(seq(1,runtime)+1500,colSums(weight.temps),col=cols(7)[i],lty=2,lwd=2)
}

image(mpa)

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
