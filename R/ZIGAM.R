# special helper functions

likelihood.fnr = function(n,det.vec,lambda,p.vec){
	det.p = prod(dbinom(det.vec,rep(n,length(det.vec)),p.vec))
	pois  = dpois(n,lambda)
	return(det.p * pois)
}

post.weight_helper = function(n,det.vec, lambda,p.vec,psi,N){
	Ns = as.matrix( min(det.vec):N )
	fns = apply(Ns,1,likelihood.fnr,det.vec,lambda,p.vec)
	Zr = psi*sum(fns) + (1-psi) * (sum(det.vec)==0)
	fn = likelihood.fnr(n,det.vec,lambda,p.vec)
	return(psi*fn/Zr)
}

log_likelihood_pois_eachsite = function(det.vec,lambda,p.vec,psi,N){
	n.site = nrow(detmat)
	nvec = max(detmat[i,]):N
	gr = apply(as.matrix(nvec),1,likelihood.fnr,det.vec=det.vec,lambda=lambda,p.vec = p.vec)
	gr = sum(gr) # given all N the probability of having data, which is actually N-mixture
	e = 1.0*(max(detmat[i,])!=0)
	logL = e * (log(psi) + log(gr))+ (1-e)*(log(1-psi + psi * gr)) # log likelihood with zero inflating 
	retrun(logL)
}

log_likelihood_pois = function(detmat,lambda,p,psi,N){
	n.site = nrow(detmat)
	logL = matrix(0,n.site,1)
	for(i in 1:n.site){
		logL[i] = log_likelihood_pois_eachsite(det.vec = detmat[i,],lambda=lambda,p.vec = p[i,],psi = psi[i],N=N) 
	}
	return(logL)
}

Hessian_sum_helper = function(detmat,lambda,p,N){
	n.site = nrow(detmat)
	lambda.sqrsumfnn = 0*detmat[,1]
	lambda.sumsqrfnn = lambda.sqrsumfnn
	psi.fnminusId0 = lambda.sqrsumfnn
	p.sumsqr=0 * detmat
	p.sum = 0 * detmat
	
	for(i in 1:n.site){
		Ns = as.matrix(max(detmat[i,]):N)
		fns = apply(Ns,1,likelihood.fnr,det.vec=detmat[i,],lambda[i],p.vec=p[i,])
		nminusmu = Ns-lambda[i]
		lambda.sqrsumfnn[i] = (sum(fns*nminusmu))^2
		lambda.sumsqrfnn[i] = sum(fns * nminusmu^2)
		psi.fnminusId0[i] = sum(fns)-(max(detmat[i,])==0)
		for(j in 1:ncol(lambda)){
			p.sumsqr[i,j] = (sum(fns * (detmat[i,j]-Ns*p[i,j])/(p[i,j]*(1-p[i,j]))))^2
			p.sum[i,j] = (1/(p[i,j]*(1-p[i,j]))^2) * sum(fns * ((detmat[i,j]-Ns*p[i,j])^2-Ns*(1-p[i,j]*p[i,j]-(1-2*p[i,j])*(detmat[i,j]-Ns * p[i,j]))))
		}
	}
	res = list(lambda.sqrsumfnn,lambda.sumsqrfnn,psi.fnminusId0,p.sumsqr,p.sum)
	return(res)
}

# helper functions for occupancy
log_likelihood_occu = function(detmat,p,psi){
	log_detli = det.vec * log(p.vec) + (1-det.vec) * log(1-det.vec)
	detli = exp(log_detli)
	Zr = psi * detli + (1-psi) * (max(detvec)==0)
	return(log(Zr))
}


occu.post.weight_helper = function(det.vec,p.vec,psi){
	log_detli = det.vec * log(p.vec) + (1-det.vec) * log(1-det.vec)
	detli = exp(log_detli)
	Zr = psi * detli + (1-psi) * (max(detvec)==0)
	return(psi * detli/Zr)
}

# These are modified from zigam(COZIGAM)

# Modified from ZIGAM.dis(COZIGAM)

## RSZIGAM.pois.R

RSZIGAM.pois <- function(formula, formula.det ,maxiter = 300, conv.crit = 1e-3,
                      size = NULL, data=list(), N,...) # data should contains detmat as det.1,det.2,det.3 etc, period should be detection period number data's formate: data$detmat should be a matrix with nrow = n.site, ncol = nperiod, data$envX, which is the second should be the environmental data at each site, the else, say data$detX.1 data$detX.2 should be the detection varible at all sites and time period 1, 2...etc., colnames should be consistent
{
  
  require(mgcv)
  gf.N.psi <- interpret.gam(formula)
  #gf.psi <- interpret.gam(formula)
  gf.det <- interpret.gam(formula.det)
  
  #y <- eval(parse(text=gf$response), envir=data)
  n.site <- nrow(data$detmat)
  period = ncol(data$detmat)
  family = poisson()
  
  d.V <- function(mu) rep.int(1, length(mu)) # the parameters as an exponential family 
  d.eta.mu <- function(mu) -1/(mu^2)
  d.f0 <- function(mu) -exp(-mu)
  den <- dpois; disp <- 1; est.disp <- FALSE
  size <- rep.int(1, n.site)
  # loglikfun <- function(y, mu, p) {
      # e <- as.numeric(y!=0)
      # sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))# likelihood with zero inflated
    # }

  
  variance <- family$variance
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  
  
  fm.psi <- as.formula(sub(gf.N.psi$response,"quasi.psi",deparse(formula)))
  fm.lambda <- as.formula(sub(gf.N.psi$response,"quasi.lambda",deparse(formula)))
  fm.p <- as.formula(sub(gf.det$response,"quasi.y",deparse(formula.det)))
  
  
  # forming data for detgams:
  detdata = matrix(nrow = 1,ncol = (1+ncol(data$envX)+ncol(data[[3]])))
  detdata = detdata[-1,]
  for(i in 1:period){ 
    y = data$detmat[,i]
    detXtemp = data[[i+2]]
	detdata = rbind(detdata, data.frame(y,data$envX,detXtemp))
  }
  ## stop here 11/14/2018 13:17 to review MATH632
  ## restart here 11/14/2018 22:20 give up MATH632
  
  lambda <- pmax(apply(data$detmat,1,mean), 5) # Poisson lambda
  psi <- rep(0.7, n.site) # occupancy psi
  # psi = runif(n.site)
  p.vec = ( 0.1*(data$detmat>=0)) # detection p initial value global 0.1
  p = matrix((p.vec),nrow = n.site*period,ncol = 1)
  # quasi.psi = matrix(runif(n.site)>0.5,nrow = n.site,ncol=1)
  quasi.lambda = pmax(apply(data$detmat,1,mean), 0.01) 
  quasi.y = data$detmat
  quasi.psi = psi
  for(i in 1:n.site){ # 
	  nvec = max(data$detmat[i,]):N
	  gr = apply(as.matrix(nvec),1,likelihood.fnr,det.vec = data$detmat[i,],lambda=lambda[i],p.vec=p.vec[i,])
      gr = sum(gr)
	  quasi.psi[i] = psi[i]*gr/(psi[i]*gr+(1-psi[i])*(sum(data$detmat[i,]!=0)==0))
  }
  norm <- 1 
  repli <- 0
  
  wg.lambda = matrix(0,n.site,1)
  wg.p = wg.lambda
  while( norm > conv.crit & repli < maxiter) { # this is the EM-PIRLS process
    
    quasi.y = data$detmat # make quasi.y the matrix form
	for(i in 1:n.site){ # again, to get quasi data of occupancy status, which is just posterior probability given all parameters
	    nvec = max(data$detmat[i,]):N
		# quasi data for occupancy status
	    gr = apply(as.matrix(nvec),1,likelihood.fnr,det.vec = data$detmat[i,],lambda=lambda[i],p.vec=p.vec[i,])
        gr = sum(gr)
	    quasi.psi[i] = psi[i]*gr/(psi[i]*gr+(1-psi[i])*(sum(data$detmat[i,]!=0)==0))
		# GAM in occupancy status has all data weight equals to 1
		
		# quasi data for Poisson lambda
		wg_tep = apply(as.matrix(nvec),1, post.weight_helper,data$detmat[i,], lambda[i],p.vec[i,],psi[i],N)
		wg.lambda[i] = sum(wg_tep)
		quasi.lambda[i] = sum(nvec * wg_tep)/wg.lambda[i]
		# weight of Poisson lambda is actually posterior probability of occupancy while pseudo data is posterior mean population size given occupy and data
		
		# quasi data for detections
		wg.p[i] = sum(nvec * wg_tep)
		quasi.y[i,] = data$detmat[i,] /quasi.lambda[i]
		# quasi y is the naive probability of having the detections given posterior Poisson mean, while weight is posterior mean of n
    }
	quasi.y = matrix(quasi.y,nrow = length(quasi.y),ncol=1) # to make quasi y a single colome
	
	
	# 11/17/2018 0:28 do not know why subscript out of bound here 11/17 12:09 Due to that formula is not well set
	G.psi <- gam(formula =  fm.psi, family = quasibinomial, fit=FALSE, data=cbind(quasi.psi, (data$envX)), ...)
	
	
	G.lambda = gam(fm.lambda,family = quasipoisson, fit = FALSE, data=data.frame(quasi.psi,data$envX),...)
	G.lambda$w = wg.lambda # change the weight in this iter, weight for the data is actually psi, see eq.9a in the technical report, here the weight is set before, see document E-step
    G.det = gam(fm.p,family = quasibinomial, fit=FALSE, data=data.frame(quasi.y,detdata[,-1]),...)
    G.det$w = rep(wg.p,period)
	fit.psi <- gam(G = G.psi) # PIRLS
    fit.lambda <- gam(G = G.lambda)
	fit.p = gam(G = G.det)
	beta.psi <- coef(fit.psi)
	beta.lambda = coef(fit.lambda)
	beta.p = coef(fit.p)
    
	lambda.old <- lambda
	p.old <- p
	psi.old = psi
	
	lambda = fit.lambda$fitted
	psi = fit.psi$fitted
	p = fit.p$fitted
	p.vec = matrix(p,nrow = n.site,ncol = period)
	
    norm <- max(abs(p-p.old), sum((lambda-lambda.old)^2),abs(psi-psi.old))
    repli <- repli + 1
    cat("iteration =", repli, "\t", "norm =", norm, "\n")
  }
 ## stop here for pokemon 11/16/2018 15:00
 ## TEST Convergence here 11/17/2018 done 12:00, results did not check, however it converges (for the simulation data set, ~270 iters)
  
  beta.psi <- coef(fit.psi)
  beta.lambda = coef(fit.lambda)
  beta.p = coef(fit.p)
  np.psi <- length(beta.psi)
  np.lambda <- length(beta.lambda)
  np.p = length(beta.p)
  sp.psi <- fit.psi$sp
  sp.lambda <- fit.lambda$sp
  sp.p = fit.p$sp
  mu.eta.val <- mu.eta(fit.lambda$linear.predictor)
  
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
  
  n.smooth <- length(G.lambda$smooth)
  Lambda.lambda <- matrix(0, np.lambda, np.lambda)
  Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
  for(k in 1:n.smooth) {
    n.S[k] <- length(G.lambda$smooth[[k]]$S)
    if(k==1) {
      Lam[[k]] <- sp.lambda[k]*G.lambda$S[[k]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp.lambda[j]*G.lambda$S[[j]]
        }
      }
    }
    else {
      Lam[[k]] <- sp.lambda[sum(n.S[1:(k-1)])+1]*G.lambda$S[[sum(n.S[1:(k-1)])+1]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp.lambda[sum(n.S[1:(k-1)])+j]*G.lambda$S[[sum(n.S[1:(k-1)])+j]]
        }
      }
    }
    first <- G.lambda$smooth[[k]]$first.para
    last <- G.lambda$smooth[[k]]$last.para
    Lambda.lambda[first:last, first:last] <- Lam[[k]]
  }

   n.smooth <- length(G.p$smooth)
   Lambda.p <- matrix(0, np.lambda, np.lambda)
   Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
   for(k in 1:n.smooth) {
     n.S[k] <- length(G.p$smooth[[k]]$S)
     if(k==1) {
       Lam[[k]] <- sp.p[k]*G.p$S[[k]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[j]*G.p$S[[j]]
         }
        }
      }
      else {
       Lam[[k]] <- sp.p[sum(n.S[1:(k-1)])+1]*G.p$S[[sum(n.S[1:(k-1)])+1]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[sum(n.S[1:(k-1)])+j]*G.p$S[[sum(n.S[1:(k-1)])+j]]
          }
        }
      }
      first <- G.p$smooth[[k]]$first.para
      last <- G.p$smooth[[k]]$last.para
      Lambda.p[first:last, first:last] <- Lam[[k]]
    }  
  
  DS.psi <- diag(eigen(Lambda.psi)$values[abs(eigen(Lambda.psi)$values)>1e-10])
  DS.lambda <- diag(eigen(Lambda.lambda)$values[abs(eigen(Lambda.lambda)$values)>1e-10])
  DS.p <- diag(eigen(Lambda.p)$values[abs(eigen(Lambda.p)$values)>1e-10])
  
  X.psi <- G.psi$X 
  X.lambda <- G.lambda$X
  X.p = G.p$X
  
  loglik <- (log_likelihood_pois(data$detmat,lambda,p.vec,psi,N)) # log-likelihood at each site, useful in calculating Hessian
  ploglik <- sum(loglik) - as.numeric(0.5*t(psi)%*%Lambda.psi%*%psi) -  as.numeric(0.5*t(lambda)%*%Lambda.lambda%*%lambda) - as.numeric(0.5*t(p)%*%Lambda.p%*%p)
  
  # stop here 13:33 11/17/2018
  # Model selection criterion
  I.theta <- matrix(0, ncol=np.psi+np.lambda+np.p, nrow=np.psi+np.lambda+np.p)  # neg Hessian at MPLE, COZIGAM has a good approximation using Laplace method to approximate the logE, including a term use this, we can derive this analytically 
  # this matrix will be block diag matrix with block to be Hessian of psi, Hessian of lambda and Hessian of p 
  # tau.lambda <- -size
  # rho.psi <- rep.int(-1, n.site)
  # rho.p = rep.int(-1,n.site*period)
  # Below is the approximation of Hessian 
  H_sum = Hessian_sum_helper(data$detmat,lambda,p.vec,N) # some useful sums, get it in single loop since I cannot avoid it.
  
  
  # Hessian block to lp of psi, no penalty yet
  Hessian.lp.psi = matrix(0,ncol = np.psi,nrow = np.psi)
  diagD = -exp(-2*loglik) * (H_sum$psi.fnminusId0)^2 * (psi*(1-psi))^2-exp(-loglik) * (H_sum$psi.fnminusId0) * (psi*(1-psi))^2 * (2*psi-1)
  D = diag(diagD)
  Hessian.lp.psi = t(X.psi) %*% D %*% X.psi
  rm(diagD)
  rm(D)
  # Hessian block to lp of lambda
  # Hessian.lp.lambda = matrix(0,ncol = np.lambda,nrow = np.lambda)
  
  diagD = -exp(-2*loglik) * psi^2 * H_sum$lambda.sqrsumfnn + exp(-loglik)*psi*H_sum$lambda.sumsqrfnn
  D = diag(diagD)
  Hessian.lp.lambda = t(X.lambda) %*% D %*% X.lambda
  rm(diagD)
  rm(D)
  
  
  # Hessian block to lp of p
  Hessian.lp.p = matrix(0,ncol = np.p,nrow = np.p)
  help_u1 = apply(H_sum$p.sumsqr,2,function(sumsqri,psi,loglik){-exp(-2*loglik) * psi^2 * sumsqri^2},psi = psi,loglik = loglik)
  help.u2 = apply(H_sum$p.sum,2,function(sumi,psi,loglik){exp(-loglik) * psi * sumi},psi = psi,loglik = loglik)
  uij = help_u1+help_u2
  rm(help_u1)
  rm(help_u2)
  diagD = matrix(uij,nrow = n.site*period,ncol = 1)
  D = diag(diagD)
  Hessian.lp.p = t(X.p)%*%D%*%X.p
  rm(diagD)
  rm(D)
  
  # add penalty 
  I.theta[1:np.psi,1:np.psi] = Hessian.lp.psi - Lambda.psi
  I.theta[(np.psi+1):(np.psi+np.lambda),(np.psi+1):(np.psi+np.lambda)] = Hessian.lp.lambda - Lambda.lambda
  I.theta[(np.psi+np.lambda+1):(np.p+np.psi+np.lambda),(np.psi+np.lambda+1):(np.p+np.psi+np.lambda)] = Hessian.lp.p - Lambda.p
  
  I.theta = -I.theta # this is the negative Hessian matrix, 
  
  logE = 0.5*determinant(DS.psi)$modulus + 0.5*determinant(DS.lambda)$modulus+ploglik + 0.5*determinant(DS.p)$modulus + 
    (np.psi+np.lambda + np.p - (nrow(DS.psi) + nrow(DS.lambda) + nrow(DS.p)))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
  attr(logE, "logarithm") <- NULL
  
  V.theta = solve(I.theta)
  V.beta.psi = V.theta[1:np.psi,1:np.psi]
  V.beta.lambda = V.theta[(np.psi+1):(np.psi+np.lambda),(np.psi+1):(np.psi+np.lambda)]
  V.beta.p = V.theta[(np.psi+np.lambda+1):(np.p+np.psi+np.lambda),(np.psi+np.lambda+1):(np.p+np.psi+np.lambda)]
  
  fit.psi$Vp = V.beta.psi
  fit.lambda$Vp = V.beta.lambda
  fit.p$Vp = V.beta.p
  
  res = list(formula = list(formula.psi=formula, formula.lambda = formula, formula.p = formula.det) ,logE=logE, ploglik=ploglik, loglik=sum(loglik),  fit.models = list(fit.psi=fit.psi, fit.lambda=fit.lambda,fit.p = fit.p),fit.values = list(psi=psi, lambda = lambda, p=p), V.beta = list( V.beta.psi=V.beta.psi, V.beta.lambda=V.beta.lambda,V.beta.p=V.beta.p), X=list(X.psi,X.lambda,X.p))
  class(res) = "RSZIGAM.Poisson.result"
  return(res)
  
}



RSZIGAM.occu <- function(formula, formula.det ,maxiter = 300, conv.crit = 1e-3,
                      size = NULL, data=list(),...) # data should contains detmat as det.1,det.2,det.3 etc, period should be detection period number data's formate: data$detmat should be a matrix with nrow = n.site, ncol = nperiod, data$envX, which is the second should be the environmental data at each site, the else, say data$detX.1 data$detX.2 should be the detection varible at all sites and time period 1, 2...etc., colnames should be consistent
{
  
  require(mgcv)
  gf.N.psi <- interpret.gam(formula)
  #gf.psi <- interpret.gam(formula)
  gf.det <- interpret.gam(formula.det)
  
  #y <- eval(parse(text=gf$response), envir=data)
  n.site <- nrow(data$detmat)
  period = ncol(data$detmat)
  family = binomial(size=1)
  
  
  fm.psi <- as.formula(sub(gf.N.psi$response,"quasi.psi",deparse(formula)))
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
  
   n.smooth <- length(G.p$smooth)
   Lambda.p <- matrix(0, np.lambda, np.lambda)
   Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
   for(k in 1:n.smooth) {
     n.S[k] <- length(G.p$smooth[[k]]$S)
     if(k==1) {
       Lam[[k]] <- sp.p[k]*G.p$S[[k]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[j]*G.p$S[[j]]
         }
        }
      }
      else {
       Lam[[k]] <- sp.p[sum(n.S[1:(k-1)])+1]*G.p$S[[sum(n.S[1:(k-1)])+1]]
       if(n.S[k]>1) {
         for(j in 2:n.S[k]) {
           Lam[[k]] <- Lam[[k]]+sp.p[sum(n.S[1:(k-1)])+j]*G.p$S[[sum(n.S[1:(k-1)])+j]]
          }
        }
      }
      first <- G.p$smooth[[k]]$first.para
      last <- G.p$smooth[[k]]$last.para
      Lambda.p[first:last, first:last] <- Lam[[k]]
    }  
  
  DS.psi <- diag(eigen(Lambda.psi)$values[abs(eigen(Lambda.psi)$values)>1e-10])
  
  DS.p <- diag(eigen(Lambda.p)$values[abs(eigen(Lambda.p)$values)>1e-10])
  
  X.psi = G.psi$X 
  
  X.p = G.p$X
  
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