#==============================================================
# Norberg et al. Eco-Evo ODE Model Functions
#    - Source this file at top of model script.
#==============================================================

#==============================================================
# Load necessary libraries
#==============================================================
library(deSolve)
library(numDeriv)

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
growth.fun<-function(rmax,TC,zi,w)
{
  rixt<-rmax*exp((-1*(TC-zi)^2)/(w^2))
  return(rixt)
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
#==============================================================
fitness.fun<-function(zi,rmax,TC,w,alphas,Nall,mort,mpa)
{
  rixt<-growth.fun(rmax=rmax,TC=TC,zi=zi,w=w)
  ints<-alphas*Nall
  if(is.null(dim(ints))) sum.ints<-sum(ints)
  else sum.ints<-apply(ints,MARGIN=2,sum)
  gixt<-rixt*(1-sum.ints)-(mort*mpa)
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
#==============================================================
dGdZ2<-function(zi,rmax,TC,w,alphas,mort,Nall,mpa)
{
  hess.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    hess.out[h]<-hessian(func=fitness.fun,x=zi[h],Nall=Nall[,h],rmax=rmax,TC=TC[h],w=w,alphas=alphas,mort=mort,mpa=mpa[h])
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
#==============================================================
dGdZ<-function(zi,rmax,TC,w,alphas,mort,Nall,mpa)
{
  grad.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    grad.out[h]<-grad(func=fitness.fun,x=zi[h],rmax=rmax,TC=TC[h],Nall=Nall[,h],w=w,alphas=alphas,mort=mort,mpa=mpa[h])
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
#==============================================================
dNdt<-function(Ni,V,zi,Di,TC,delx,rmax,w,alphas,mort,Nall,mpa)
{
  popdy<-Ni*fitness.fun(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,mort=mort,Nall=Nall,mpa=mpa) # Population dynamics component
  genload<-.5*V*Ni*dGdZ2(zi=zi,rmax=rmax,TC=TC,w=w,Nall=Nall,alphas=alphas,mort=mort,mpa=mpa) # Genetic load component
  dispersal<-Di*dNdx2(Ni=Ni,delx=delx) # Dispersal component
  popchange<-popdy+genload+dispersal
  popchange[popchange<(-1*Ni+10^-6)]<-(-1*Ni[popchange<(-1*Ni+10^-6)]+10^-6)  # Checking if population density falls below Nmin
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
#==============================================================  
dZdt<-function(Ni,V,zi,Di,TC,delx,rmax,w,alphas,mort,Nall,mpa,Nmin)
{               
  q<-qfun(Ni=Ni)
  directselect<-q*V*dGdZ(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,
                               Nall=Nall,mort=mort,mpa=mpa)                     # Directional selection component
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
  
  return(list("dZbardt"=dZbar.dt,"Ecology"=eco,"Evolution"=evo,"props"=p_all))
  
}
