# vertual data formation
#setwd("C:/Users/yshen99/Documents/GitHub/RSZIGAM/R")
source("RSZIGAM.pois.R")
source('misc.R')
meshgrid = function (xrange, yrange){
  ncol = length(xrange)
  nrow = length(yrange)
  X.x = matrix(nrow = nrow,ncol = ncol)
  X.y = matrix(nrow = nrow,ncol = ncol)
  for (i in 1:ncol){
    X.y[,i] = yrange
  }
  for (i in 1:nrow){
    X.x[i,] = xrange
  }
  
  X = list(X.x,X.y)
  return(X)
}
require(e1071)

#setwd('E:/UW Lab jobs/Community Assembly/ANN Poisson/vertual data')
n.site = 400
n.period = 20
set.seed(12345)
env.1 = 5*runif(n.site)
env.2 = 5*runif(n.site)
env.3 = 5*runif(n.site)
det.1 = 5*runif(n.period)
P.xy = meshgrid(det.1,env.1)

#lambda = 4* dnorm(env.1,mean = 0.5,sd = 0.3) * dnorm(env.2,mean = 0.75,sd = 0.2)
#lambda = dnorm(env.1,mean = 0.4,sd=0.15)+dnorm(env.1,mean=0.7,sd = 0.1) + dnorm(env.2,mean = 0.75,sd=0.2)+ dnorm(env.3,mean = 0.75,sd=0.2)
lambda = env.1+env.2+env.3
#p = sigmoid(as.matrix( P.xy[[1]])) * 
#  dnorm(as.matrix( P.xy[[2]]))/max(as.matrix( dnorm(P.xy[[2]])))
#p = sigmoid(as.matrix( P.xy[[1]]) + 
#  dnorm(as.matrix( P.xy[[2]]))/max(as.matrix( dnorm(P.xy[[2]]))))
p = sigmoid(0.4*P.xy[[1]]-0.3*P.xy[[2]])
#psi = sigmoid(env.3*10-5)*dnorm(env.1,sd=0.7)/(max(dnorm(env.1,sd=0.7)))
#psi = sigmoid(env.3+dnorm(env.1,0.3,sd=0.2)/(max(dnorm(env.1,0.3,sd=0.2))))
psi = sigmoid(0.8*env.1+0.3*env.2-env.3)

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

vdata = data.frame(lambda = lambda,psi = psi,Ni,occu.r,p,detection)
parameters = data.frame(lambda = lambda,psi = psi,p)
realpopu = data.frame(N = Ni,occu = !occu.r)
env = data.frame(env.1,env.2,env.3)

# write.csv(parameters,'parameters.csv',row.names = F)
# write.csv(detection,'countingdata.csv',row.names = F)
# write.csv(realpopu,'population.csv',row.names = F)
# write.csv(env,'env.csv',row.names = F)
# write.csv(P.xy[[1]],'detection_coef.csv',row.names = F)

detX = t( matrix(det.1,nrow = 20,ncol = 30))
data.test = list(detmat = detection,envX=env)
for(i in 1:n.period + 2){
  data.test[[i]] = data.frame(det.1 = rep(det.1[i-2],n.site))
}


RES.simu=RSZIGAM.pois(y~s(env.1,bs="cr",k=-1)+s(env.2,bs="cr",k=-1)+s(env.3,bs="cr",k=-1),y~s(env.1,bs="cr",k=-1)+s(det.1,bs="cr",k=-1),data=data.test,N=100,maxiter = 100)
# no error yet, debug mode

RES.simu.linear=RSZIGAM.pois(y~env.1+env.2+env.3,y~env.1+det.1,data=data.test,N=100,maxiter = 100)


envnew = data.frame(env.1=env.1,env.2=.2,env.3=.3)

plot(env.1,predict(RES.simu$fit.models$fit.lambda,envnew))
