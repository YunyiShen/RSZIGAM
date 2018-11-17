# special helper functions

likelihood.fnr = function(n,det.vec,lambda,p.vec){
	det.p = prod(dbinom(det.vec,rep(n,length(det.vec)),p.vec))
	pois  = dpois(n,lambda)
	return(det.p * pois)
}

log_likelihood = function(detmat,lambda,p,psi,N){
	n.site = nrow(detmat)
	logL = 0
	for(i in 1:n.site){
		nvec = max(detmat[i,]):N
		gr = apply(as.matrix(nvec),1,likelihood.fnr,det.vec = detmat[i,],lambda=lambda[i],p.vec = p[i,])
		gr = sum(gr) # given all N the probability of having data, which is actually N-mixture
		e = 1.0*(max(detmat[i,])!=0)
		logL = logL + e * (log(psi[i]) + log(gr))+ (1-e)*(log(1-psi[i] + psi[i] * gr)) # log likelihood with zero inflating 
	}
	return(logL)
}

post.weight_helper = function(n,det.vec, lambda,p.vec,psi,N){
	Ns = as.matrix( min(det.vec):N )
	fns = apply(Ns,1,likelihood.fnr,det.vec,lambda,p.vec)
	Zr = psi*sum(fns) + (1-psi) * (sum(det.vec)==0)
	fn = likelihood.fnr(n,det.vec,lambda,p.vec)
	return(psi*fn/Zr)
}

# These are modified from zigam(COZIGAM)

# Modified from ZIGAM.dis(COZIGAM)

## RSZIGAM.dis.R

RSZIGAM.dis <- function(formula, formula.det ,maxiter = 20, conv.crit = 1e-3,
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
  loglikfun <- function(y, mu, p) {
      e <- as.numeric(y!=0)
      sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))# likelihood with zero inflated
    }

  
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
  # psi <- rep(0.7, n.site) # occupancy psi
  psi = runif(n.site)
  p.vec = ( 0.1*(data$detmat>=0)) # detction p
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
  while( norm > conv.crit & repli < maxiter) { # this is the EM-PIRLS process hopefully 
    
  quasi.y = data$detmat
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
  
  loglik <- log_likelihood(data$detmat,lambda,p.vec,psi,N) # log-likelihood
  ploglik <- loglik - as.numeric(0.5*t(b.psi)%*%Lambda.psi%*%b.psi) -  as.numeric(0.5*t(b.lambda)%*%Lambda.lambda%*%b.lambda) - as.numeric(0.5*t(b.p)%*%Lambda.p%*%b.p)
  
  # stop here 13:33 11/17/2018
  # Model selection criterion
  I.theta <- matrix(0, ncol=np1+np2, nrow=np1+np2)
  tau.mu <- -size
  rho.p <- rep.int(-1, n)
  a <- 1-p+p*den(0, mu)
  tau.mu[y==0] <- ( -p*mu.eta.val/disp/variance(mu)/a*(den(0,mu)-size*mu*den(0,mu)*
                                                         (variance(mu)*d.eta.mu(mu)+d.V(mu)/mu.eta.val)*mu.eta.val/variance(mu)+size*(1-p)*mu*d.f0(mu)/a) )[y==0]
  rho.p[y==0] <- ( (1-den(0,mu))/(a^2)*(-a*(1-2*p)-(1-den(0,mu))*p*(1-p)) )[y==0]
  
  G.tau.mu <- diag(as.vector(mu.eta.val*tau.mu))
  G.rho.p <- diag(as.vector(p*(1-p)*rho.p))
  
  I.theta[1:np1,1:np1] <- t(X1) %*% G.tau.mu %*% X1 - Lambda1
  I.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)] <- t(X2) %*% G.rho.p %*% X2 - Lambda.lambda
  I.theta <- -I.theta
  
  logE <- 0.5*determinant(DS1)$modulus + 0.5*determinant(DS2)$modulus+ploglik +
    (np1+np2-(nrow(DS1)+nrow(DS2)))/2*log(2*pi)-0.5*determinant(I.theta)$modulus
  attr(logE, "logarithm") <- NULL
  
  V.theta <- solve(I.theta)
  V.beta <- V.theta[1:np1, 1:np1]
  V.gamma <- V.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)]
  
  fit.gam$Vp <- V.beta; fit.lr$Vp <- V.gamma
  
  res <- list(formula=formula, logE=logE, ploglik=ploglik, loglik=loglik, fit.gam=fit.gam, fit.lr=fit.lr,
              mu=mu, p=p, psi=psi, dispersion=disp, V.beta=V.beta, V.gamma=V.gamma, X1=X1, X2=X2, family=family)
  res
  
}

