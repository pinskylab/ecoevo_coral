source("Norberg_Functions.r")

nsp<-4
size<-50
times<-seq(0,80,by=1)
temps<-20+.1*seq(1,size,length.out=size)
Topt1<-temps
Topt2<-temps
Topt3<-temps
Topt4<-temps
Nmin<-10^-6
maxw<-7

# sppstate<-matrix(c(rep(.01,5),rep(Nmin,2*size-10),rep(.01,5)),nrow=nsp,ncol=size,byrow=T)
sppstate<-matrix(c(rep(.2,size*nsp)),nrow=nsp,ncol=size,byrow=T)
traitstate<-matrix(c(rep(Topt1,size),rep(Topt2,size),rep(Topt3,size),rep(Topt4,size)),nrow=nsp,ncol=size,byrow=T)
# 
# traitstate<-matrix(c(rep(Topt1,size),rep(Topt2,size),rep(Topt3,size),rep(Topt4,size)),nrow=nsp,ncol=size,byrow=T)

sptype=c(1,1,1,2)
mpa<-matrix(rep(1,size*nsp),nrow=nsp,ncol=size,byrow=T)
mpa[sptype==1,c(5:10,20:30,40:45)]<-.75
mpa[sptype==2,c(20:30)]<-1.25


allnames<-c("Comp","Stress","Comp_opt1","Stress_opt2","temps")

allstate<-c(sppstate[1,],sppstate[2,],sppstate[3,],sppstate[4,],traitstate[1,],traitstate[2,]
            ,traitstate[3,],traitstate[4,],temps)
#temps<-rep(23,size)
# parms<-list(
#   
#   V=c(.3,.3,.3,.3),
#   D=c(0.1,0.1,.1,.1),
#   nsp=nsp,
#   rmax=c(1,.3,.5,1),
#   alphas=matrix(1,nrow=nsp,ncol=nsp),
#   m=c(.15,.1,.2,.3),
#   w=maxw*c(1,1,1,1),
#   Nmin=10^-6,
#   deltax=1
# 
# )

coral_trait_stoch<-function(t,y, parms,size,nsp){
  
  ts<-matrix(NA,nrow=(2*nsp*size+size),ncol=t)
  dspout<-ts
  ts[,1]<-y
  
  for(k in 2:t)
  {
    print(k)
    V=c(.5,.5,.5,.5)
    D=c(0.1,0.1,.1,.1)
    nsp=nsp
    rmax=c(.85,.2,.75,1)
    alphas=matrix(1,nrow=nsp,ncol=nsp)
    m=c(.15,.05,.17,.25)
    w=maxw*c(1,1,1,1)
    
    Nmin=10^-6
    deltax=1
    
    spps<-matrix(ts[(1:(size*nsp)),k-1],nrow=nsp,ncol=size,byrow=T)
    traits<-matrix(ts[(size*nsp+1):(length(y)-size),k-1],nrow=nsp,ncol=size,byrow=T)
    temps<-ts[(length(y)-(size-1)):length(y),k-1]
    
    pcat<-rbinom(1,1,.15)
    
    m<-(1-pcat)*m+pcat*c(.75,.1,.55,.3)
    
    dspp<-spps
    dtraits<-traits
    
    for(i in 1:nsp)
    { 
      dspp[i,]<-dNdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,])
      dtraits[i,]<-dZdt(Ni=spps[i,],zi=traits[i,],Nall=spps,rmax=rmax[i],V=V[i],D=D[i],w=w[i],mort=m[i],TC=temps,alphas=alphas[i,],delx=1,mpa=mpa[i,],Nmin=Nmin)
    }
    
    dtemps<-rep(0.00,size)#+rnorm(size,0,.1) 
    
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
    dspout[,k]<-dsp
    ts[,k]<-ts[,k-1]+dsp
    
  }
  
  return(list(ts=ts,dsp=dspout))

}

out<-coral_trait_stoch(t=500,y=allstate, parms=parms,size=size,nsp=4)

imageCols<-colorRampPalette(c("blue","red","gold"))
layout(matrix(c(1,2,3,4,
                5,6,7,8,
                9,9,10,11),nrow=3,ncol=4,byrow=T))
for(i in 1:9)
{
  if(i<5){
  image(x=seq(1,500),y=seq(1:size),z=t(out$ts[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
        ylab="Location on Reef",zlim=c(0,.8))
  }
  else{
    image(x=seq(1,500),y=seq(1:size),z=t(out$ts[((i-1)*size+1):(i*size),]),col=imageCols(24),xlab="Time",
          ylab="Location on Reef",zlim=c(20,25))
  }
}

plot(seq(1,500),(colSums(out$ts[((1-1)*size+1):(1*size),])/50),ylim=c(0,1),type="l",
     col="darkgreen",yaxs="i",bty="l",ylab="Hard Coral Cover",xlab="Time")
lines(seq(1,500),(colSums(out$ts[((2-1)*size+1):(2*size),])/50),col="red")
lines(seq(1,500),(colSums(out$ts[((3-1)*size+1):(3*size),])/50),col="blue")
lines(seq(1,500),(colSums(out$ts[((1-1)*size+1):(3*size),])/50),col="black",lwd=2)

plot(seq(1,500),(colSums(out$ts[((4-1)*size+1):(4*size),])/50),ylim=c(0,1),type="l",
     col="darkgreen",yaxs="i",bty="l",ylab="MacroAlgae Cover",xlab="Time",lwd=2)


