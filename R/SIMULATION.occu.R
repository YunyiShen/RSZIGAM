source("RSZIGAM.occu.R")
source("misc.R")

n.site = 900
n.period = 10
set.seed(42)
env.1 = 5*runif(n.site)
env.2 = 5*runif(n.site)
env = data.frame(env.1,env.2)

det.1 = matrix( 5*runif(n.period*n.site),nrow = n.site,ncol=n.period )


psi = 3 * (dnorm(env.1,mean = 1,sd=.5) + dnorm(env.1,mean=3,sd = 1)) + 
         3*dnorm(env.2,mean = 3,sd=1)
psi = sigmoid(psi)
z = runif(n.site)<=psi
p = sigmoid(apply(det.1,2,function(det,env){1*det-1*env},env=env.1))

det = matrix(runif(n.period*n.site),n.site,n.period)<=p
detection = apply(det,2,function(det,z){det*z},z=z)

data.test = list(detmat = detection,envX=env)
for(i in 1:n.period + 2){
  data.test[[i]] = data.frame(det.1 = det.1[,i-2])
}


RES.simu=RSZIGAM.occu(y~s(env.1,bs="tp",k=-1)+s(env.2,bs="tp",k=-1),
                      y~s(env.1,bs="tp",k=-1)+s(env.2,bs="tp",k=-1)+s(det.1,bs="tp",k=-1),
                      data=data.test,maxiter = 50)

#check norm difference
norm(as.matrix(RES.simu$fit.values$psi-psi),"1")/norm(as.matrix(psi),"1")
norm(as.matrix(RES.simu$fit.values$p-p),"1")/norm(as.matrix(p),"1")

# check reaction curves
# psi
envnew = data.frame(env.1=env.1,env.2=.2)
plot(env.1,sigmoid(predict(RES.simu$fit.models$fit.psi,envnew)))
plot(env.1,sigmoid(3 * (dnorm(env.1,mean = 1,sd=.5) + dnorm(env.1,mean=3,sd = 1))))

envnew = data.frame(env.1=.1,env.2=env.2)
plot(env.2,(predict(RES.simu$fit.models$fit.psi,envnew)))
plot(env.2,2*dnorm(env.2,mean = 3,sd=1))

#p
detenvnew = data.frame(env.1=.1,env.2=.2,det.1=5*runif(1000))
plot(detenvnew$det.1,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$det.1,detenvnew$det.1)

detenvnew = data.frame(env.1=env.1,env.2=.2,det.1=0*runif(100))
plot(detenvnew$env.1,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$env.1,-detenvnew$env.1)

detenvnew = data.frame(env.1=.1,env.2=env.2,det.1=.2)
plot(detenvnew$env.2,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$env.2,0*detenvnew$det.1)

raster::plot(raster::raster(matrix(RES.simu$fit.values$psi,30,30)))
raster::plot(raster::raster(matrix(psi,30,30)))
raster::plot(raster::raster(matrix((psi-RES.simu$fit.values$psi)/psi,30,30)),zlim=c(-1,1))

