# vertual data formation
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
n.site = 30
n.period = 20
set.seed(0)
env.1 = runif(n.site)
env.2 = runif(n.site)
env.3 = runif(n.site)
det.1 = runif(n.period)
P.xy = meshgrid(det.1,env.1)

#lambda = 4* dnorm(env.1,mean = 0.5,sd = 0.3) * dnorm(env.2,mean = 0.75,sd = 0.2)
lambda = dnorm(env.1,mean = 0.4,sd=0.15)+dnorm(env.1,mean=0.7,sd = 0.1) + dnorm(env.2,mean = 0.75,sd=0.2)
#p = sigmoid(as.matrix( P.xy[[1]])) * 
#  dnorm(as.matrix( P.xy[[2]]))/max(as.matrix( dnorm(P.xy[[2]])))
p = sigmoid(as.matrix( P.xy[[1]]) + 
  dnorm(as.matrix( P.xy[[2]]))/max(as.matrix( dnorm(P.xy[[2]]))))
#psi = sigmoid(env.3*10-5)*dnorm(env.1,sd=0.7)/(max(dnorm(env.1,sd=0.7)))
psi = sigmoid(env.3+dnorm(env.1,0.3,sd=0.2)/(max(dnorm(env.1,0.3,sd=0.2))))

Ni = apply(as.matrix(lambda),1,rpois,n=1)
detection = matrix(nrow =n.site, ncol = n.period)
occu.p = runif(n.site)
occu.r = occu.p>psi
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


RSZIGAM.pois(y~s(env.1,bs="cr",k=4)+s(env.2,bs='cr',k=3)+s(env.3,bs='cr',k=3),y~s(env.1,bs='cr',k=3)+s(env.2,bs='cr',k=3)+s(env.3,bs='cr',k=3)+s(det.1,bs='cr',k=3),data=data.test,N=100,maxiter = 1)
# no error yet, debug mode
