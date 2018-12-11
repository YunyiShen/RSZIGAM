source("RSZIGAM.pois.R")
source('misc.R')

require(e1071)

logit = function(p){
  linearp = log(p/1-p)
  return(linearp)
}


n.site = 400
n.period = 20
set.seed(12345)
env.1 = 5*runif(n.site)
env.2 = 5*runif(n.site)

det.1 = matrix( 5*runif(n.period*n.site),nrow = n.site,ncol=n.period )


lambda = 3 * (dnorm(env.1,mean = 1,sd=.5) + dnorm(env.1,mean=3,sd = 1) + 
         2*dnorm(env.2,mean = 3,sd=1.5))

psi = sigmoid(dnorm(env.2,2,sd=1)/max(dnorm(env.2,2,sd=1))+
                dnorm(env.1,3,sd=1)/(max(dnorm(env.1,3,sd=1))))

p = sigmoid(apply(det.1,2,function(det,env){1*det-1*env},env=env.1))

Ni = apply(as.matrix(lambda),1,rpois,n=1)
detection = matrix(nrow =n.site, ncol = n.period)
occu.p = runif(n.site)
occu.r = occu.p>psi # not occu
for (i in 1:n.site){
  if(occu.r[i]){
    detection[i,] = rep(0,n.period)
    next
  }
  for (j in 1:n.period){
    detection[i,j] = rbinom(1,size = Ni[i],prob = p[i,j])
  }
}

env = data.frame(env.1,env.2)

data.test = list(detmat = detection,envX=env)
for(i in 1:n.period + 2){
  data.test[[i]] = data.frame(det.1 = det.1[,i-2])
}


RES.simu=RSZIGAM.pois(y~s(env.1,bs="cr",k=-1)+s(env.2,bs="cr",k=-1),
                      y~s(env.1,bs="cr",k=-1)+s(env.2,bs="cr",k=-1)+s(det.1,bs="cr",k=-1),
                      data=data.test,N=100,maxiter = 50)

#check norm difference
norm(as.matrix(RES.simu$fit.values$psi-psi),"1")/norm(as.matrix(psi),"1")
norm(as.matrix(RES.simu$fit.values$lambda-lambda),"1")/norm(as.matrix(lambda),"1")
norm(as.matrix(RES.simu$fit.values$p-p),"1")/norm(as.matrix(p),"1")

# check reaction curves
# lambda
envnew = data.frame(env.1=env.1,env.2=.2)
plot(env.1,exp(predict(RES.simu$fit.models$fit.lambda,envnew)))
plot(env.1,3 * (dnorm(env.1,mean = 1,sd=.5) + dnorm(env.1,mean=3,sd = 1)))

envnew = data.frame(env.1=.1,env.2=env.2)
plot(env.2,exp(predict(RES.simu$fit.models$fit.lambda,envnew)))
plot(env.2,2*dnorm(env.2,mean = 3,sd=1.5))

# psi
envnew = data.frame(env.1=(env.1),env.2=.2)
plot(env.1,(predict(RES.simu$fit.models$fit.psi,envnew)))
plot(env.1,dnorm(env.1,3,sd=1)/(max(dnorm(env.1,3,sd=1))))

envnew = data.frame(env.1=.1,env.2=env.2)
plot(env.2,predict(RES.simu$fit.models$fit.psi,envnew))
plot(env.2,dnorm(env.2,2,sd=1)/(max(dnorm(env.1,2,sd=1))))


detenvnew = data.frame(env.1=.1,env.2=.2,det.1=5*runif(100))
plot(detenvnew$det.1,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$det.1,detenvnew$det.1)

detenvnew = data.frame(env.1=env.1,env.2=.2,det.1=0*runif(100))
plot(detenvnew$env.1,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$env.1,-detenvnew$env.1)

detenvnew = data.frame(env.1=.1,env.2=env.2,det.1=.2)
plot(detenvnew$env.2,predict(RES.simu$fit.models$fit.p,detenvnew))
plot(detenvnew$env.2,0*detenvnew$det.1)
