source("Norberg_Functions.r")

nsp<-2
size<-50
times<-seq(0,500,by=1)
Topt1<-20
Topt2<-25
Nmin=10^-6

sppstate<-matrix(c(rep(.01,5),rep(Nmin,90),rep(.01,5)),nrow=nsp,ncol=size,byrow=T)
traitstate<-matrix(c(rep(Topt1,size),rep(Topt2,size)),nrow=2,ncol=size,byrow=T)

allnames<-c("spp1","spp2","opt1","opt2","temps")


temps<-20+.1*seq(1,50,length.out=size)

allstate<-c(sppstate[1,],sppstate[2,],traitstate[1,],traitstate[2,],temps)

parms<-list(
  
  V=.5,
  D=0.1,
  nsp=nsp,
  rmax=1,
  alphas=matrix(1,nrow=nsp,ncol=nsp),
  m=.1,
  w=5,
  
  Nmin=10^-6,
  deltax=1
  
  )




coral_trait<-function(t,y, parms,size,nsp,names){
  

  
  with(as.list(parms), {
  
    spps<-matrix(y[1:(size*nsp)],nrow=nsp,ncol=size,byrow=T)
    traits<-matrix(y[(size*nsp+1):(length(y)-size)],nrow=nsp,ncol=size,byrow=T)
    temps<-y[(length(y)-(size-1)):length(y)]

    dspp<-spps
    dtraits<-traits
    for(i in 1:nsp)
    { 
      dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax,V=V,D=D,w=w,mort=m,TC=temps,alphas=alphas[i,],delx=1)
      dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax,V=V,D=D,w=w,mort=m,TC=temps,alphas=alphas[i,],delx=1)
    }
    
    dtemps<-rep(0.02,size)
    
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
    return(list(dsp))
  })

  
}

coral_trait(t=times,y=allstate,parms=parms,size=size,nsp=nsp)

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


