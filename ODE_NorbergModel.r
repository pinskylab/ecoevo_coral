source("Norberg_Functions.r")

nsp<-2                        # how many species in model?
size<-20                      # how many reefs in model?
times<-seq(0,200,by=1)        # How many time steps?
mid<-25                       # mean temperature across all reefs at start of simulation.
range<-5                      # range of temperatures across reefs at start of simulation

temps<-generate.temps(size=size,mid=mid,range=range,temp.scenario="linear")
sppstate<-generate.state(size=size,nsp=nsp,dens=0.1)
traitstate<-generate.traits(nsp=nsp,size=size,mid=mid,range=range,trait.scenario="u.const")
allstate<-c(as.vector(t(sppstate)),as.vector(t(traitstate)),temps)

allnames<-c("spp1","spp2","opt1","opt2","temps")

parms<-list(
  
  V=.5,
  D=0.1,
  nsp=nsp,
  rmax=1,
  alphas=matrix(1,nrow=nsp,ncol=nsp),
  m=.1,
  w=5,
  mpa=.9,
  annual.temp.change<-.02,
  maxtemp<-30,
  
  Nmin=10^-6,
  deltax=1
  
  )

coral_trait<-function(t,y, parms,size,nsp,names,temp.change=c("const","linear","sigmoid")){
    
  with(as.list(parms), {
  
    spps<-matrix(y[1:(size*nsp)],nrow=nsp,ncol=size,byrow=T)
    traits<-matrix(y[(size*nsp+1):(length(y)-size)],nrow=nsp,ncol=size,byrow=T)
    temps<-y[(length(y)-(size-1)):length(y)]

    dspp<-spps
    dtraits<-traits
    for(i in 1:nsp)
    { 
      dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax,V=V,D=D,w=w,mort=m,
                     TC=temps,alphas=alphas[i,],delx=1,mpa=mpa)
      dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax,V=V,D=D,w=w,mort=m,
                        TC=temps,alphas=alphas[i,],delx=1,mpa=mpa)
    }
    
    dtemps<-switch(temp.change,
                   const=rep(0,size),
                   linear=rep(annual.temp.change,size),
                   sigmoid=rep((annual.temp.change*mean(temps)*(1-(mean(temps)/maxtemp))),size))
    
    dsp<-NA
    j<-1
    while(j <=nsp)
    {
      dsp<-c(dsp,dspp[j,])
      j<-j+1
    }
    dsp<-dsp[-1]
    j<-1
    while(j <=nsp)
    {
      dsp<-c(dsp,dtraits[j,])
      j<-j+1
    }
    dsp<-c(dsp,dtemps)
   # dsp<-c(as.vector(t(dspp)),as.vector(t(dtraits)),dtemps)
    return(list(dsp))
  })

  
}

coral_trait(t=times,y=allstate,parms=parms,size=size,nsp=nsp,temp.change="constant")

out <- ode.1D(y = allstate, t = times, func = coral_trait, parms = parms,
              nspec = (nsp*2)+1, names = allnames, dimens=size,size=size,nsp=nsp)

delx<-1
Distance <- seq(from = 0.5, by = delx, length.out = size)

#plot individual species
# par(mfrow=c(2,2))
windows()

image(out, which=c("spp1","spp2","opt1","opt2","temps"), grid = Distance,
      xlab = "time, days", ylab = "Distance on Reef, km",
      main = c("Coral 1 density By Reef Location","Coral 2 density By Reef Location",
               "Coral 1 Trait By Reef Location","Coral 2 Trait By Reef Location","Temperature"),legend=T)
# #windows()
# image(out, which="spp2", grid = Distance,
#       xlab = "time, days", ylab = "Distance on Reef, km",
#       main = "Coral 2 density By Reef Location",legend=T)
# 
# #windows()
# image(out, which="opt1", grid = Distance,
#       xlab = "time, days", ylab = "Distance on Reef, km",
#       main = "Coral 1 Trait By Reef Location",legend=TRUE)
# #windows()
# image(out, which="opt2", grid = Distance,
#       xlab = "time, days", ylab = "Distance on Reef, km",
#       main = "Coral 2 Trait By Reef Location",legend=TRUE)


