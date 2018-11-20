
# These are modified from zigam(COZIGAM)

# Modified from ZIGAM.dis(COZIGAM)

# data should contains detmat as det.1,det.2,det.3 etc, period should be detection period number data's formate: data$detmat should be a matrix with nrow = n.site, ncol = nperiod, data$envX, which is the second should be the environmental data at each site, the else, say data$detX.1 data$detX.2 should be the detection varible at all sites and time period 1, 2...etc., colnames should be consistent

RSZIGAM.occu <- function(formula, formula.det ,maxiter = 300, conv.crit = 1e-3,
                      size = NULL, data=list(),...)
{
  source("misc.R")
  require(mgcv)
  #gf.N.psi <- interpret.gam(formula)
  gf.psi <- interpret.gam(formula)
  gf.det <- interpret.gam(formula.det)
  
  #y <- eval(parse(text=gf$response), envir=data)
  n.site <- nrow(data$detmat)
  period = ncol(data$detmat)
  family = binomial(size=1)
  
  
  fm.psi <- as.formula(sub(gf.psi$response,"quasi.psi",deparse(formula)))
  fm.p <- as.formula(sub(gf.det$response,"y",deparse(formula.det)))
  
  
  # forming data for detgams:
  detdata = matrix(nrow = 1,ncol = (1+ncol(data$envX)+ncol(data[[3]])))
  detdata = detdata[-1,]
  for(i in 1:period){ 
    y = data$detmat[,i]
    detXtemp = data[[i+2]]
	detdata = rbind(detdata, data.frame(y,data$envX,detXtemp))
  }
  
  
  lambda <- pmax(apply(data$detmat,1,mean), 5) # Poisson lambda
  psi <- rep(0.7, n.site) # occupancy psi
  # psi = runif(n.site)
  p.vec = ( 0.1*(data$detmat>=0)) # detection p initial value global 0.1
  p = matrix((p.vec),nrow = n.site*period,ncol = 1)
  # quasi.psi = matrix(runif(n.site)>0.5,nrow = n.site,ncol=1)
  quasi.psi = psi
  norm <- 1 
  repli <- 0
  
  wg.lambda = matrix(0,n.site,1)
  wg.p = wg.lambda
  while( norm > conv.crit & repli < maxiter) { # this is the EM-PIRLS process
    
    quasi.y = data$detmat # make quasi.y the matrix form
	for(i in 1:n.site){ # again, to get quasi data of occupancy status, which is just posterior probability given all parameters
		# quasi data for occupancy status
	    quasi.psi[i] = occu.post.weight_helper(det.vec=data$detmat[i,],p.vec=p.vec[i,],psi)
		# GAM in occupancy status has all data weight equals to 1
    }
	wg.p[i] = quasi.psi
	y = matrix(detmat,nrow = length(detmat),ncol=1) # to make quasi y a single colome
	
	
	
	G.psi <- gam(formula =  fm.psi, family = quasibinomial, fit=FALSE, data=cbind(quasi.psi, (data$envX)), ...)
	
	
	# change the weight in this iter, weight for the data is actually psi, see eq.9a in the technical report, here the weight is set before, see document E-step
    G.det = gam(fm.p,family = quasibinomial, fit=FALSE, data=detdata,...)
    G.det$w = rep(wg.p,period)
	fit.psi <- gam(G = G.psi) # PIRLS
	fit.p = gam(G = G.det)
	beta.psi = coef(fit.psi)
	beta.p = coef(fit.p)
    
	p.old <- p
	psi.old = psi
	
	psi = fit.psi$fitted
	p = fit.p$fitted
	p.vec = matrix(p,nrow = n.site,ncol = period)
	
    norm <- max(abs(p-p.old), sum((lambda-lambda.old)^2),abs(psi-psi.old))
    repli <- repli + 1
    cat("iteration =", repli, "\t", "norm =", norm, "\n")
  }
   
  beta.psi <- coef(fit.psi)
  
  beta.p = coef(fit.p)
  np.psi <- length(beta.psi)
  
  np.p = length(beta.p)
  sp.psi <- fit.psi$sp
  
  sp.p = fit.p$sp
  
  
  ## penalty for three GAMs
  n.smooth <- length(G.psi$smooth)
  Lambda.psi <- matrix(0, np.psi, np.psi)
  Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
  for(k in 1:n.smooth) {
    n.S[k] <- length(G.psi$smooth[[k]]$S)
    if(k==1) {
      Lam[[k]] <- sp.psi[k]*G.psi$S[[k]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp.psi[j]*G.psi$S[[j]]
        }
      }
    }
    else {
      Lam[[k]] <- sp.psi[sum(n.S[1:(k-1)])+1]*G.psi$S[[sum(n.S[1:(k-1)])+1]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp.psi[sum(n.S[1:(k-1)])+j]*G.psi$S[[sum(n.S[1:(k-1)])+j]]
        }
      }
    }
    first <- G.psi$smooth[[k]]$first.para
    last <- G.psi$smooth[[k]]$last.para
    Lambda1[first:last, first:last] <- Lam[[k]]
  }
  
   n.smooth <- length(G.det$smooth)
   Lambda.p <- matrix(0, np.lambda, np.lambda)
   Lam <- list()
   n.S <- numeric(n.smooth) # penalty matrix
   for(k in 1:n.smooth) {
     n.S[k] <- length(G.det$smooth[[k]]$S)
     if(k==1) {
       Lam[[k]] <- sp.p[k]*G.det$S[[k]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[j]*G.det$S[[j]]
         }
        }
      }
      else {
       Lam[[k]] <- sp.p[sum(n.S[1:(k-1)])+1]*G.det$S[[sum(n.S[1:(k-1)])+1]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[sum(n.S[1:(k-1)])+j]*G.det$S[[sum(n.S[1:(k-1)])+j]]
          }
        }
      }
      first <- G.det$smooth[[k]]$first.para
      last <- G.det$smooth[[k]]$last.para
      Lambda.p[first:last, first:last] <- Lam[[k]]
    }  
  
  DS.psi <- diag(eigen(Lambda.psi)$values[abs(eigen(Lambda.psi)$values)>1e-10])
  
  DS.p <- diag(eigen(Lambda.p)$values[abs(eigen(Lambda.p)$values)>1e-10])
  
  X.psi = G.psi$X 
  
  X.p = G.det$X
  
  loglik <- log_likelihood_occu(data$detmat,p.vec,psi) # log-likelihood
  ploglik <- loglik - as.numeric(0.5*t(psi)%*%Lambda.psi%*%psi) - as.numeric(0.5*t(p)%*%Lambda.p%*%p)
  
  # stop here 13:33 11/17/2018
  # Model selection criterion
  I.theta <- matrix(0, ncol=np.psi+np.p, nrow=np.psi+np.p)  # neg Hessian at MPLE, COZIGAM has a good approximation using Laplace method to approximate the logE, including a term use this, we can derive this analytically 
  # this matrix will be block diag matrix with block to be Hessian of psi, Hessian of lambda and Hessian of p 
  rho.psi <- rep.int(-1, n.site)
  rho.p = rep.int(-1,n.site*period)
  # Below is the approximation of Hessian 
  # Hessian block to lp of psi, no penalty yet
  Hessian.lp.psi = matrix(0,ncol = np.psi,nrow = np.psi)
  # Hessian block to lp of p
  Hessian.lp.p = matrix(0,ncol = np.p,nrow = np.p)
  
  # add penalty 
  I.theta[1:np.psi,1:np.psi] = Hessian.lp.psi - Lambda.psi
  I.theta[(np.psi+1):(np.p+np.psi),(np.psi+1):(np.p+np.psi)] = Hessian.lp.p - Lambda.p
  
  I.theta = -I.theta # this is the negative Hessian matrix, 
  
  logE = 0.5*determinant(DS.psi)$modulus +ploglik + 0.5*determinant(DS.p)$modulus + 
    (np.psi+np.lambda + np.p - (nrow(DS.psi) + nrow(DS.p)))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
  attr(logE, "logarithm") <- NULL
  
  V.theta = solve(I.theta)
  V.beta.psi = V.theta[1:np.psi,1:np.psi]
  V.beta.p = V.theta[(np.psi+1):(np.p+np.psi),(np.psi+1):(np.p+np.psi)]
  
  fit.psi$Vp = V.beta.psi
  fit.p$Vp = V.beta.p
  
  res = list(formula = list(formula.psi=formula, formula.p = formula.det) ,logE=logE, ploglik=ploglik, loglik=loglik,  fit.models = list(fit.psi=fit.psi,fit.p = fit.p),fit.values = list(psi=psi, p=p), V.beta = list( V.beta.psi=V.beta.psi,V.beta.p=V.beta.p), X=list(X.psi,X.p))
  class(res) = "RSZIGAM.Occupancy.result"
  return(res)
  
}