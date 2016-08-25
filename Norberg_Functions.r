library(deSolve)
library(numDeriv)


growth.fun<-function(rmax,TC,zi,w)
{
  rixt<-rmax*exp((-1*(TC-zi)^2)/(w^2))
  return(rixt)
}

fitness.fun<-function(zi,rmax,TC,w,alphas,Nall,mort,mpa)
{
  rixt<-growth.fun(rmax=rmax,TC=TC,zi=zi,w=w)
  ints<-alphas*Nall
  if(is.null(dim(ints))) sum.ints<-sum(ints)
  else sum.ints<-apply(ints,MARGIN=2,sum)
  gixt<-rixt*(1-sum.ints)-(mort*mpa)
  return(gixt)
}

dNdx<-function(Ni,delx) # can use for both dNdx and dzdx
{
  dNdx<-diff(c(Ni,0))/delx
}

dNdx2<-function(Ni,delx) # can use for both dNdx and dzdx
{
  d1<-dNdx(Ni=c(0,Ni),delx=delx)
  d2<-diff(d1)/delx
  return(d2 )
}

dZdx<-function(zi,delx) # can use for both dNdx and dzdx
{
  dZdx<-diff(c(zi,zi[length(zi)]))/delx
}

dZdx2<-function(zi,delx) # can use for both dNdx and dzdx
{
  d1<-dZdx(zi=c(zi[1],zi),delx=delx)
  d2<-diff(d1)/delx
  return(d2 )
}


my.hessian<-function(zi,rmax,TC,w,alphas,N,mort,Nall,mpa)
{
  hess.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    hess.out[h]<-hessian(func=fitness.fun,x=zi[h],Nall=Nall[,h],rmax=rmax,TC=TC[h],w=w,alphas=alphas,mort=mort,mpa=mpa[h])
  }
  return(hess.out)
}

my.grad<-function(zi,rmax,TC,w,alphas,N,mort,Nall,mpa)
{
  grad.out<-rep(0,length(zi))
  for(h in 1:length(zi))
  {
    grad.out[h]<-grad(func=fitness.fun,x=zi[h],rmax=rmax,TC=TC[h],Nall=Nall[,h],w=w,alphas=alphas,mort=mort,mpa=mpa[h])
  }
  return(grad.out)
}



dNdt<-function(Ni,V,zi,Di,TC,delx,rmax,w,alphas,mort,Nall,mpa)
{
  popdy<-Ni*fitness.fun(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,mort=mort,Nall=Nall,mpa=mpa)
  genload<-.5*V*Ni*my.hessian(zi=zi,rmax=rmax,TC=TC,w=w,Nall=Nall,alphas=alphas,N=Ni,mort=mort,mpa=mpa)
  dispersal<-Di*dNdx2(Ni=Ni,delx=delx)
  popchange<-popdy+genload+dispersal
  popchange[popchange<(-1*Ni+10^-6)]<-(-1*Ni[popchange<(-1*Ni+10^-6)]+10^-6)
  return(popchange)
}

qfun<-function(Ni,Nmin=10^-6)
{
  qixt<-max(0,
            1-(Nmin/(max(Nmin,2*Ni))))
  return(qixt)
}

  
dZdt<-function(V,Ni,zi,TC,rmax,w,alphas,mort,D,delx,Nall,mpa,Nmin)
{
  q<-qfun(Ni=Ni)
  directselect<-q*V*my.grad(zi=zi,rmax=rmax,TC=TC,w=w,alphas=alphas,
                               Nall=Nall,mort=mort,mpa=mpa)
  geneflow<-D*(dZdx2(zi=zi,delx=delx)+2*dNdx(Ni=log(max(Ni,Nmin)),
                                             delx=delx)*dZdx(zi=zi,delx=delx))
  traitchange<-directselect+geneflow
  return(traitchange)
  
}

#==============================================
# Example code to test functions
#==============================================
# 
# temps<-20+.1*(seq(1,30,by=1))
# nsp<-2
# Topt<-20
# traits<-rep(Topt,(nsp*30))
# V<-1
# D<-0
# delx<-.1
# w<-5
# rmax<-1
# alphas<-matrix(c(1,1,1,1),byrow=T,nrow=2,ncol=2)
# m<-.1


# hessian(func=fitness.fun,x=20,rmax=rmax,TC=22.5,w=w,alphas=alphas[1,],N=.1,mort=m)
# grad(func=fitness.fun,x=20,rmax=rmax,TC=22.5,w=w,alphas=alphas[1,],N=.1,mort=m)
# # 
# my.hessian(zi=traits[1:30],rmax=rmax,TC=temps,w=w,alphas=alphas[1,],
#            Nall=rbind(rep(.1,30),rep(.1,30)),mort=m)
# my.grad(zi=traits[1:30],rmax=rmax,TC=temps,w=w,alphas=alphas[1,],N=rep(.1,30),mort=m)
# 
# gixt<-fitness.fun(zi=traits[1:30],rmax=rmax,TC=rep(21.5,30),w=w,alphas=alphas[1,],
#             N=rep(.1,30),m=m)



# 
# dNdt(Ni=rep(.01,30),V=V,zi=temps,Di=D,TC=temps,delx=delx,
#      rmax=1,alphas=alphas[1,],mort=m,w=w)


