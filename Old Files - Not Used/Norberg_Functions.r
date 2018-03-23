#==============================================================
# Norberg et al. Eco-Evo ODE Model Functions
#    - Source this file at top of model script.
#==============================================================

# ==============================================================
# Load necessary libraries
# ==============================================================
# install.packages(c("deSolve","numDeriv","doSNOW","doParallel","foreach",
#                    "rbenchmark","fitdistrplus","logitnorm","matrixcalc",
#                    "mvtnorm"))
library(deSolve)
library(numDeriv)
library(doSNOW)
library(doParallel)
library(foreach)
library(rbenchmark)
library(fitdistrplus)
library(logitnorm)
library(matrixcalc)
library(mvtnorm)
#library(mail)

#==============================================================
# Growth rate function (Equation (3) from supplemental material
#   to Norberg et al. 2012). Determines growth rate as a 
#   function of species trait (optimum temperature) and local
#   environmental condition (reef temperature).
# Parameters:
#   rmax: Maximum growth rate for species
#   TC:   Current temperature of reef
#   zi:   Current optimum temperature for species
#   w:    Temperature tolerance
#==============================================================
growth.fun<-function(rmax,TC,zi,w,spp="C",growth="normal")
{

  rixt<-(rmax/sqrt(2*w^2*pi))*exp((-1*(TC-zi)^2)/(2*w^2))
  if(growth=="Huey"){
    CTmax<-zi+w
    rixt[zi>=TC]<-(rmax/sqrt(2*w^2*pi))*(exp((-1*(TC-zi)^2)/(2*w)))[zi>=TC]
    rixt[zi<TC]<-(rmax/sqrt(2*w^2*pi))*(1-((TC-zi)/(zi-CTmax))^2)[zi<TC]
    for(i in 1:length(rixt)) rixt[i]<-max(0,rixt[i])
  }
  if(growth=="skew"){
  rixt<-skew.norm(rmax,TC,zi,w)
  }
  if(spp=="MA") rixt<-rep(.49,length(zi))*rmax
  return(rixt)
}

#==============================================================
# Mortality rate function proposed as an alternative to
#   constant mortality rate from Norberg et al. Determines
#   mortality rate as a function of species trait (optimum 
#   temperature) and local environmental condition (reef
#   temperature).
# Parameters:
#   TC:   Current temperature of reef
#   zi:   Current optimum temperature for species
#   w:    Temperature tolerance
#==============================================================
mortality.fun<-function(TC,zi,w,rmax=NA,spp="C",mpa,growth="normal",almort)
{
  
  mixt<-1-exp((-1*(TC-zi)^2)/(w^2))
  mixt[zi>=TC]<-0
  mixt[mixt<0.03]<-0.03
  if(growth=="Huey"& spp!="MA"){
    CTmax<-zi+w
    mixt[zi>=TC]<-(rmax/sqrt(2*w^2*pi)-(rmax/sqrt(2*w^2*pi))*(exp((-1*(TC-zi)^2)/(2*w))))[zi>=TC]
    mixt[zi<TC]<-(1-(1-(rmax/sqrt(2*w^2*pi))*((TC-zi)/(zi-CTmax))^2))[zi<TC]
    for(i in 1:length(mixt)) mixt[i]<-min(1,mixt[i])
  }
  if(spp=="MA"){
    # mixt<-rep(0,length(TC))
    mixt[mpa==1] <- 0.3
    mixt[mpa!=1] <- almort[mpa!=1]

    #return(mixt)
  }
  #print("Coral")
  #print(mixt)
  return(mixt)
}



#==============================================================
# Fitness function (Equation (2) from supplemental material
#   to Norberg et al. 2012). Determines fitness as a 
#   function of local growth rate, mortality rate, and species
#   interactions (competition for space).
# Calls growth.fun()
# Parameters:
#   rmax:   Maximum growth rate for species
#   TC:     Current temperature of reef
#   zi:     Current optimum temperature for species
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   Nall:   Vector of abundances for all species
#   mort:   Mortality rate
#   mpa:    Effect of MPA on mortality rate of species
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================
fitness.fun<-function(zi,rmax,TC,w,alphas,Nall,mort,mpa,
                      mortality.model=c("constant","tempvary"),spp,growth="normal")
{
  rixt<-growth.fun(rmax=rmax,TC=TC,zi=zi,w=w,spp=spp,growth=growth)
  if(mortality.model=="tempvary"){

    mort<-mortality.fun(TC=TC,zi=zi,rmax=rmax,w=w,spp=spp,mpa=mpa,
                        growth=growth,almort=mort)
    
  } 

  ints<-alphas*Nall
  if(is.null(dim(ints))) sum.ints<-sum(ints)
  else sum.ints<-apply(ints,MARGIN=2,sum)
  gixt<-rixt*(1-sum.ints)-(mort)
  return(gixt)
}

#==============================================================
# Function to calculate the partial derivitive of species
#   density over space. This function is used in the last
#   expression of equation (1a) and the last expression of
#   equation (1b) from Norberg et al. 2012 supplemental
#   material. Note that when used in equation (1b), need to 
#   enter log(Ni) for Ni. For use in equation (1a), it is 
#   called within dNdx2(), described below.
# Parameters:
#   Ni:   Vector of densities of species i along reef 
#   delx: Interval across which to calculate derivitive  
#==============================================================
dNdx<-function(Ni,delx)
{
  dNdx<-diff(c(Ni,0))/delx
}

#==============================================================
# Function to calculate the second partial derivitive of
#   species density over space. This function is used in the
#   last expression of equation (1a) from Norberg et al. 2012
#   supplemental material.
# Calls dNdx()
# Parameters:
#   Ni:   Vector of densities of species i along reef 
#   delx: Interval across which to calculate derivitive  
#==============================================================
dNdx2<-function(Ni,delx) # can use for both dNdx and dzdx
{
  d1<-dNdx(Ni=c(0,Ni),delx=delx)
  d2<-diff(d1)/delx
  return(d2 )
}

#==============================================================
# Function to calculate the partial derivitive of trait
#   value over space. This function is used in the last 
#   expression of equation (1b) from Norberg et al. 2012 
#   supplemental material. For use in equation (1b), it is 
#   called within dNdx2(), described below, as well as called
#   alone.
# Parameters:
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   delx: Interval across which to calculate derivitive  
#==============================================================
dZdx<-function(zi,delx) # can use for both dNdx and dzdx
{
  dZdx<-diff(c(zi,zi[length(zi)]))/delx
}

#==============================================================
# Function to calculate the second partial derivitive of trait
#   value over space. This function is used in the last 
#   expression of equation (1b) from Norberg et al. 2012 
#   supplemental material. 
# Calls dZdx()
# Parameters:
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   delx: Interval across which to calculate derivitive  
#==============================================================
dZdx2<-function(zi,delx)
{
  d1<-dZdx(zi=c(zi[1],zi),delx=delx)
  d2<-diff(d1)/delx
  return(d2 )
}

#==============================================================
# Function to calculate the second partial derivitive of growth
#   rate across changes in trait space. This function is used 
#   the genetic load component of equation (1a) in the 
#   Norberg et al. 2012 supplemental material. 
# Calls fitness.fun(), growth.fun()
# Parameters:
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   rmax:   Maximum growth rate for species
#   TC:     Current temperature of reef
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   Nall:   Vector of abundances for all species
#   mort:   Mortality rate
#   mpa:    Effect of MPA on mortality rate of species
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================
dGdZ2<-function(zi,rmax,TC,w,alphas,mort,Nall,mpa,mortality.model,spp,growth="normal")
{
  hess.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    hess.out[h]<-hessian(func=fitness.fun,x=zi[h],Nall=Nall[,h],
                         rmax=rmax,TC=TC[h],w=w,alphas=alphas,
                         mort=mort[h],mpa=mpa[h],mortality.model=mortality.model,spp=spp,growth=growth)
  }
  return(hess.out)
}

#==============================================================
# Function to calculate the partial derivitive of growth
#   rate across changes in trait space. This function is used 
#   the directionsal selection component of equation (1b) in 
#   the Norberg et al. 2012 supplemental material. 
# Calls fitness.fun(), growth.fun()
# Parameters:
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   rmax:   Maximum growth rate for species
#   TC:     Current temperature of reef
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   Nall:   Vector of abundances for all species
#   mort:   Mortality rate
#   mpa:    Effect of MPA on mortality rate of species
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================
dGdZ<-function(zi,rmax,TC,w,alphas,mort,Nall,mpa,mortality.model,spp,growth="normal")
{
  grad.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    grad.out[h]<-grad(func=fitness.fun,x=zi[h],rmax=rmax,
                      TC=TC[h],Nall=Nall[,h],w=w,alphas=alphas,
                      mort=mort[h],mpa=mpa[h],
                      mortality.model=mortality.model,method="simple",spp=spp,growth=growth)
  }
  return(grad.out)
}

#==============================================================
# Function to calculate the partial derivitive of density
#   across time. Equation (1a) in the Norberg et al. 2012 
#   supplemental material. 
# Calls fitness.fun(), growth.fun(), dNdx(), dNdx2(), dGdZ2(),
#   dGdZ()
# Parameters:
#   Ni:   vector of current densities foe species i
#   V:    Genetic variation
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   Di:     Dispersal rate for species i
#   TC:     Current temperature of reef
#   delx: Interval across which to calculate derivitive for
#           spatial derivitives
#   rmax:   Maximum growth rate for species
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   mort:   Mortality rate
#   Nall:   Vector of abundances for all species
#   mpa:    Effect of MPA on mortality rate of species
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================
dNdt<-function(Ni,V,zi,Di,TC,delx,rmax,w,alphas,mort,Nall,mpa,
               mortality.model,spp,growth="normal")
{  
  popdy<-Ni*fitness.fun(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,
                        mort=mort,Nall=Nall,mpa=mpa,
                        mortality.model=mortality.model,spp=spp,growth=growth) # Population dynamics component
  genload<-.5*V*Ni*dGdZ2(zi=zi,rmax=rmax,TC=TC,w=w,Nall=Nall,
                         alphas=alphas,mort=mort,mpa=mpa,
                         mortality.model=mortality.model,spp=spp,growth=growth) # Genetic load component
  dispersal<-Di*dNdx2(Ni=Ni,delx=delx) # Dispersal component
  popchange<-popdy+genload+dispersal
  popchange[(Ni+popchange)<10^-6 | is.na(popchange)]<- -Ni[(Ni+popchange)<10^-6 | is.na(popchange)]+10^-6
  #popchange[popchange<(-1*Ni+10^-6)]<-(-1*Ni[popchange<(-1*Ni+10^-6)]+10^-6)  # Checking if population density falls below Nmin
  return(popchange)
}

#==============================================================
# Function to calculate the partial derivitive of density
#   across time. Equation (1a) in the Norberg et al. 2012 
#   supplemental material. 
# Calls fitness.fun(), growth.fun(), dNdx(), dNdx2(), dGdZ2(),
#   dGdZ()
# Parameters:
#   Ni:   vector of current densities foe species i
#   V:    Genetic variation
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   Di:     Dispersal rate for species i
#   TC:     Current temperature of reef
#   delx: Interval across which to calculate derivitive for
#           spatial derivitives
#   rmax:   Maximum growth rate for species
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   mort:   Mortality rate
#   Nall:   Vector of abundances for all species
#   mpa:    Effect of MPA on mortality rate of species
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================
dNdt.vec<-function(data.in,V,Di,delx,rmax,w,alphas,mort,mpa,
                    mortality.model,spp,growth="normal")
{
  
  Ni<-data.in[,1]
  zi<-data.in[,2]
  TC<-data.in[,3]
  Nall<-data.in[,4:6]
  
  popdy<-Ni*fitness.fun(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,
                        mort=mort,Nall=Nall,mpa=mpa,
                        mortality.model=mortality.model,spp=spp,growth=growth) # Population dynamics component
  genload<-.5*V*Ni*dGdZ2(zi=zi,rmax=rmax,TC=TC,w=w,Nall=Nall,
                         alphas=alphas,mort=mort,mpa=mpa,
                         mortality.model=mortality.model,spp=spp,growth=growth) # Genetic load component
  dispersal<-Di*dNdx2(Ni=Ni,delx=delx) # Dispersal component
  popchange<-popdy+genload+dispersal
  popchange[(Ni+popchange)<10^-6 | is.na(popchange)]<- -Ni[(Ni+popchange)<10^-6 | is.na(popchange)]+10^-6
  #popchange[popchange<(-1*Ni+10^-6)]<-(-1*Ni[popchange<(-1*Ni+10^-6)]+10^-6)  # Checking if population density falls below Nmin
  return(popchange)
}


#==============================================================
# Function to prevent directional selection of virtually
#   extinct populations and enhance numerical stability. 
#   Described in supplemental material of Norberg et al (2012).
# Parameters:
#   Ni:   vector of current densities foe species i
#   Nmin: Minimum value for density at any given location
#==============================================================
qfun<-function(Ni,Nmin=10^-6)
{
  qixt<-max(0,
            1-(Nmin/(max(Nmin,2*Ni))))
  return(qixt)
}

#==============================================================
# Function to calculate the partial derivitive of trait value
#   across time. Equation (1b) in the Norberg et al. 2012 
#   supplemental material. 
# Calls fitness.fun(), growth.fun(), dNdx(), dNdx2(), dGdZ2(),
#   dGdZ()
# Parameters:
#   Ni:   vector of current densities foe species i
#   V:    Genetic variation
#   zi:   Vector of optimal temperatures for species i along 
#             reef 
#   Di:     Dispersal rate for species i
#   TC:     Current temperature of reef
#   delx: Interval across which to calculate derivitive for
#           spatial derivitives
#   rmax:   Maximum growth rate for species
#   w:      Temperature tolerance
#   alphas: Species interaction matrix
#   mort:   Mortality rate
#   Nall:   Vector of abundances for all species
#   mpa:    Effect of MPA on mortality rate of species
#   Nmin:   Minimum density allowed (required for log values)
#   mortality.model  Should constant or temperature varying
#                     mortality be used?
#==============================================================  
dZdt<-function(Ni,V,zi,Di,TC,delx,rmax,w,alphas,mort,Nall,mpa,Nmin,
               mortality.model,spp,growth="normal")
{ 
  
  q<-qfun(Ni=Ni)
  directselect<-q*V*dGdZ(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,
                               Nall=Nall,mort=mort,mpa=mpa,
                         mortality.model=mortality.model,spp=spp,growth=growth)                     # Directional selection component
  geneflow<-Di*(dZdx2(zi=zi,delx=delx)+2*dNdx(Ni=log(max(Ni,Nmin)),
                                             delx=delx)*dZdx(zi=zi,delx=delx))  # Gene flow component
  traitchange<-directselect+geneflow
  return(traitchange)
  
}

#==============================================================
# Function to partition the ecological and evolutionary
#   processes responsible for local trait changes. From the
#   supplemental material of Norberg et al. (2012). Calculates
#   values at every point in space and time.
# Parameters:
#   Zall:  Vector of traits for all species at site j across
#             time. Rows for species, columns for time.
#   Nall:  Vector of densities for all species at site j
#             across time. Rows for species, columns for time.
#==============================================================
dZbardt<-function(Zall,Nall)
{
  
  # Ecology component
  p_all<-t(t(Nall)/colSums(Nall))             # proportional densities for each species across time
  p_all.add<-cbind(p_all[,1],p_all)           # add a column at start of time series so that derivitive at
                                              #   time 0 is 0.
  dpdt<-t(apply(p_all.add,MARGIN=1,FUN=diff))
  eco<-colSums(Zall*dpdt)
  
  # Evolution component
  Zall.add<-cbind(Zall[,1],Zall)
  dZ.dt<-t(apply(Zall.add,MARGIN=1,FUN=diff))
  evo<-colSums(p_all*dZ.dt)
  
  dZbar.dt<-eco+evo
  
  return(list("dZbardt"=dZbar.dt,"Ecology"=abs(eco),"Evolution"=abs(evo),"props"=p_all))
  
}

#======================================================================
# The following functions are not from the Norberg et al 2012 paper,
#     but are used to generate starting conditions for species density,
#     optimal traits, and temperatures across the reef.
#======================================================================

#======================================================================
# Function to generate temperatures across the reef. Temperatures can
#     take on several patterns, including uniform, linear increase, or
#     randomized.
# Parameters:
#     size:  Size of the reef (number of cells)
#     mid:   Mean of the temperature range on reef
#     range: Range of temperatures across reef
#     temp.scenario: Temperature scenario. Uniform has same (mid value)
#                       temperature at all locations. Linear has a 
#                       linearly increasing temperature centered on
#                       mid value. Random draws a random temperature
#                       for each cell from a normal distribution
#                       with mean equal to mid value, and sd equal to
#                       half of range value.
#======================================================================
generate.temps<-function(size,mid=25,range=2.5,temp.scenario=c("uniform","linear","random"))
{
  temps<-switch(temp.scenario,
                uniform=rep(mid,size),
                linear=seq((mid-(range)),(mid+range),length.out=size),
                random=seq((mid-(range)),(mid+range),length.out=size)[sample(size=size,x=seq(1,size),replace=F)])
  
  return(temps)
}

#======================================================================
# Function to generate traits across the reef. Traits can
#     take on several patterns, including uniform, linear increase, or
#     randomized. Traits represent the optimum temperature for species
#     i at each location.
# Parameters:
#     nsp:   Number of species in model
#     size:  Size of the reef (number of cells)
#     mid:   Mean of the temperature range on reef
#     range: Range of temperatures across reef
#     temps: Temperatures at each reef location
#     trait.scenario: Trait scenario. u.const has unique traits for
#                       each species, which are constant across reef. 
#                       perfect.adapt has perectly adapted organisms,
#                       whose traits are equal to the temperature on
#                       they experience on the reef. Same.constant 
#                       produces the same trait for all species at all
#                       locations.
#======================================================================
generate.traits<-function(nsp,size,mid,range,temps,trait.scenario=c("u.const","perfect.adapt",
                                                                    "same.constant"))
{
  if(trait.scenario%in%c("u.const","perfect.adapt",
                         "same.constant")){
    if(trait.scenario=="u.const")
    {
      trts<-seq((mid-(range/4)),(mid+range/4),length.out=(nsp+2))[-c(1,nsp+2)]
      traits<-matrix(trts,nrow=nsp,ncol=size,byrow=F)
    }
    
    if(trait.scenario=="same.const")
    {
      traits<-matrix(mid,nrow=nsp,ncol=size,byrow=F)
    }
    
    if(trait.scenario=="perfect.adapt")
    {
      traits<-matrix(temps,nrow=nsp,ncol=size,byrow=T)
    }
    
  }
  
  else stop("Invalid scenario provided")
  
  return(traits)
}

#======================================================================
# Function to generate starting densities for all species at each
#     location on the reef. Currently produces constant densities
#     across species and space.
# Parameters:
#     nsp:   Number of species in model
#     size:  Size of the reef (number of cells)
#     dens:  Starting density
#======================================================================
generate.state<-function(size,nsp,dens=.01,random=F)
{
  state<-matrix(dens,nrow=nsp,ncol=size)
  if(random){
    state<-matrix(runif((nsp*size),10E-6,.33),nrow=nsp,ncol=size)
  }
  return(state)
}

#=====================================================================
# Function to calculate Shannon diversity 
#
#=====================================================================
shannon.div<-function(N)
{
  sdiv<-0
  for(i in 1:length(N))
  {
    sdiv<-sdiv+(N[i]*log(N[i]))
  }
  sdiv<--1*sdiv
  return(sdiv)
}

#====================================================================
# Gauss Error Function
#====================================================================
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


#===================================================================
# Skew- Normal growth
#===================================================================
skew.norm<-function(rmax,TC,zi,w,lambda=-2.7)
{
  growth<-rmax*exp(-1*(TC-zi)^2/w^2)*(1+erf(lambda*(TC-zi)/w))
  
  return(growth)
  
}


#===================================================================
# Function to set MPA design relative to protection strategy
#===================================================================
setMPA<-function(temps,Nall,sptype,strategy=c("hot","cold","highcoral","lowcoral",
                                              "portfolio","random","none"),size,amount=.2,
                 priordata=NULL,monitortime=20)
{
  mpa<-rep(0,size)
  if(amount==0) return(mpa)
  corals<-colSums(Nall[sptype==1,])
  ncoral<-sum(sptype==1)
  if(strategy=="hot")
  {
    ind<-order(temps,decreasing=TRUE)[1:(amount*size)]
  }
  if(strategy=="cold")
  {
    ind<-order(temps)[1:(amount*size)]
  }
  if(strategy=="hotcold")
  {
    ind<-order(temps)[c((1:(amount*size/2)),(size:((size-amount*size/2+1))))]
    
  }
  if(strategy=="space")
  {
   ind<-round(seq(1,size,length.out=amount*size)) 
    
  }
  if(strategy=="highcoral")
  {
    ind<-order(corals,decreasing=TRUE)[1:(amount*size)]
  }
  if(strategy=="lowcoral")
  {
    ind<-order(corals)[1:(amount*size)]
  }
  if(strategy=="random")
  {
    ind<-sample(1:size,amount*size)
  }
  if(strategy=="portfolio")
  {
    ind<-eff.portfolio(Nall=priordata,sptype=sptype,size=size,amount=amount,
                       monitortime = monitortime)
    
  }
  if(strategy=="portfolioGreedy")
  {
    ind<-greedyport(Nall=priordata,sptype=sptype,size=size,amount=amount,
                   monitortime = monitortime)
  }
  if(strategy!="none") mpa[ind]<-1
  return(mpa)
}


#=================================================================
# Function to select an efficient portfolio
#=================================================================

eff.portfolio<-function(Nall,sptype=c(1,1,2),
                        size,amount=.2,monitortime,samples=100000)
{
  ncoral<-sum(sptype==1)
  maxtime<-ncol(Nall)
  if(ncoral==2){
    dat1<-Nall[(1:size),(maxtime-(monitortime-1)):maxtime]
    dat2<-Nall[(size+1):(2*size),(maxtime-(monitortime-1)):maxtime]
    dat<-dat1+dat2
  }
  if(ncoral!=2) stop("Need to add functionality to eff.portfolio beyond 2 coral species")
  
  ExR<-colMeans(t(dat))
  VcovR<-cov(t(dat))
  RtoSBest<-0
  RtoSCurrent<-0
  portfolios<-matrix(NA,nrow=2,ncol=amount*size)
  for(i in 1:samples)
  {
    w<-rep(0,size)
    portfolio<-sample(1:size,amount*size,replace=F)
    #print(portfolio)
    w[portfolio]<-1/(amount*size)
    
    Sig<-t(w)%*%VcovR%*%w
    
    RtoSCurrent<-(w%*%ExR)/Sig
    portfolios[2,]<-portfolio
    if(RtoSCurrent>RtoSBest|i==1)
    {
      portfolios[1,]<-portfolio
      RtoSBest<-RtoSCurrent
    }
  }
  
  # best<-which(RtoS==max(RtoS))
  #print(best)
  return(portfolios[1,])  
  

}


greedyport<-function(Nall,sptype=c(1,1,2),
                    size,amount=.2,monitortime)
{
  
  ncoral<-sum(sptype==1)
  maxtime<-ncol(Nall)
  if(ncoral==2){
    dat1<-Nall[(1:size),(maxtime-(monitortime-1)):maxtime]
    dat2<-Nall[(size+1):(2*size),(maxtime-(monitortime-1)):maxtime]
    dat<-dat1+dat2
  }
  if(ncoral!=2) stop("Need to add functionality to eff.portfolio beyond 2 coral species")
  
  ExR<-colMeans(t(dat))
  VcovR<-cov(t(dat))
  RtoSBest<-0
  RtoSCurrent<-0
  portfolioCurr<-0
  portfolioBest<-0
  
  pcombs<-combn(size,2)
  for(j in 1:ncol(pcombs))
  {
    portfolio<-pcombs[,j]
    w<-rep(0,size)
    #portfolio<-sample(1:size,amount*size,replace=F)
    #print(portfolio)
    w[portfolio]<-1/(amount*size)
    
    Sig<-t(w)%*%VcovR%*%w
    
    RtoSCurrent<-(w%*%ExR)/Sig
    portfolioCurr<-portfolio
    if(RtoSCurrent>RtoSBest|j==1)
    {
      portfolioBest<-portfolioCurr
      RtoSBest<-RtoSCurrent
    }
  }
  k<-2
  
  while(k<(amount*size))
  {
    portfolioLast<-portfolioBest
    RtoSLast<-RtoSBest
    remainingsites<-(1:size)[-c(portfolioBest)]
    for(z in remainingsites)
    {
      portfolio<-c(portfolioLast,z)

      w<-rep(0,size)
      #portfolio<-sample(1:size,amount*size,replace=F)
      #print(portfolio)
      w[portfolio]<-1/(amount*size)
      #print(w)
      Sig<-t(w)%*%VcovR%*%w
      #print(Sig)
      RtoSCurrent<-(w%*%ExR)/Sig
      #print(RtoSCurrent)
      portfolioCurr<-portfolio
      if(RtoSCurrent>RtoSBest|z==remainingsites[1])
      {
        portfolioBest<-portfolioCurr
        RtoSBest<-RtoSCurrent
        
      }
      
    }

    k<-k+1
  }
  

  
  return(portfolioBest)
  
  
}



#====================================================================
# Function to denote extinctions across entirety of reef
# Parameters:
#   Nall: matrix of species densities across space and time
#   nsp: number of species
#   threshold: Extinction threshold
#====================================================================
extinctrisk<-function(Nall,nsp=3,threshold=.001,size)
{
  sp.extinct<-rep(NA,nsp)
  for(i in 1:(nsp))
  {
    print(length(Nall[i,]))
    spdat<-sum(Nall[i,]<threshold)
    sp.extinct[i]<-spdat==size
  }
  return(sp.extinct)
}

extinctrisk2<-function(Nall,nsp=3,threshold=.001)
{
  sp.extinct<-apply(Nall,MARGIN=1,FUN=mean)<threshold
  return(sp.extinct)
}




