# special helper functions

liklihood.fnr = function(n,det.vec,lambda,p.vec){
	det.p = prod(dbinom(det.vec,rep(n,length(det.vec)),p.vec))
	pois  = dpois(n,lambda)
	return(det.p * pois)
}

post.weight_helper = function(n,det.vec, lambda,p.vec,psi,N){
	Ns = as.matrix( min(det.vec):N )
	fns = apply(Ns,2,likelihood.fnr,det.vec,lambda,p.vec)
	Zr = psi*sum(fns) + (1-psi) * (sum(det.vec)==0)
	fn = likelihood.fnr(n,det.vec,lambda,p.vec)
	return(psi*fn/Zr)
}

# These are modified from zigam(COZIGAM)

zigam <- function(formula, maxiter = 20, conv.crit = 1e-3,
                  size = NULL, log.tran = FALSE, family, data=list(), ...)
{
  
  if (is.character(family))
    fam <- eval(parse(text = family))
  if (is.function(family))
    fam <- family()
  if (fam$family == "gaussian" | fam$family == "Gamma") {
    zigam.res <- ZIGAM.cts(formula, log.tran = log.tran, family = fam, data=data, ...)
    attr(zigam.res, "family.type") <- "continuous"
  }
  else if (fam$family == "poisson" | fam$family == "binomial") {
    zigam.res <- ZIGAM.dis(formula, maxiter, conv.crit, size = size, family = fam, data=data, ...)
    attr(zigam.res, "family.type") <- "discrete"
  }
  else stop("family not recognized")
  attr(zigam.res, "constraint") <- "none"
  
  invisible(zigam.res)
  
}



# ZIGAM.dis(COZIGAM)

## ZIGAM.dis.R

ZIGAM.dis <- function(formula, formula.det ,maxiter = 20, conv.crit = 1e-3,
                      size = NULL, data=list(), N,...) # data should contains detmat as det.1,det.2,det.3 etc, period should be detection period number data's formate: data$detmat should be a matrix with nrow = nsite, ncol = nperiod, data$envX, which is the second should be the environmental data at each site, the else, say data$detX.1 data$detX.2 should be the detection varible at all sites and time period 1, 2...etc., colnames should be consistent
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
  size <- rep.int(1, n)
  loglikfun <- function(y, mu, p) {
      e <- as.numeric(y!=0)
      sum((1-e)*log(1-p+p*dpois(y,mu,log=FALSE))+e*(log(p)+dpois(y,mu,log=TRUE)))# likelihood with zero inflated
    }

  
  variance <- family$variance
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  
  
  fm.psi <- as.formula(sub(gf.N.psi$response,"psi",deparse(formula)))
  fm.n <- as.formula(sub(gf.N.psi$response,"n",deparse(formula)))
  fm.p <- as.formula(sub(gf.det$response,"y",deparse(formula)))
  
  # forming data for detgams:
  detdata = list()
  for(i in 1:period){
    y = data$detmat[,i]
    detXtemp = data[[i+2]]
	detdata[[i]]=data.frame(y,data$envX,detXtemp)
  }
  ## stop here 11/14/2018 13:17 for review MATH632
  
  
  mu <- pmax(y, 0.01)
  p <- rep(0.7, n)
  psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0))
  norm <- 1; repli <- 0
  while( norm > conv.crit & repli < maxiter) { # this is the EM-PIRLS process hopefully 
    
    psi <- p*den(y, mu)/(p*den(y, mu)+(1-p)*(y==0)) # this is E step, meanwhile, it calculated the weight for PIRLS
	# seems we need another E step here regarding the Latent N, and P, using different weight to deal with it:
	for(n in 0:N){
	n_pois = n * size
	
	
	
	
	
	}
	
	# then M step
    G1 <- gam(fm1, family=family, fit=FALSE, data=data, ...)
    G2 <- gam(fm2, family=quasibinomial, fit=FALSE, data=data, ...)
    G1$w <- psi*size # change the weight in this iter, weight for the data is actually psi, see eq.9a in the technical report
    fit.gam <- gam(G = G1) # seems this is the PIRLS work, done by gam 
    fit.lr <- gam(G = G2)
    b <- coef(fit.gam)
    g <- coef(fit.lr)
    
    mu.old <- mu; p.old <- p
    mu <- fit.gam$fitted
    p <- fit.lr$fitted
    norm <- max(abs(p-p.old), sum((mu-mu.old)^2))
    repli <- repli + 1
    cat("iteration =", repli, "\t", "norm =", norm, "\n")
  }
  
  b1 <- fit.gam$coef; b2 <- fit.lr$coef
  np1 <- length(b1); np2 <- length(b2)
  sp1 <- fit.gam$sp; sp2 <- fit.lr$sp
  mu.eta.val <- mu.eta(fit.gam$linear.predictor)
  
  n.smooth <- length(G1$smooth)
  Lambda1 <- matrix(0, np1, np1)
  Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
  for(k in 1:n.smooth) {
    n.S[k] <- length(G1$smooth[[k]]$S)
    if(k==1) {
      Lam[[k]] <- sp1[k]*G1$S[[k]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp1[j]*G1$S[[j]]
        }
      }
    }
    else {
      Lam[[k]] <- sp1[sum(n.S[1:(k-1)])+1]*G1$S[[sum(n.S[1:(k-1)])+1]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp1[sum(n.S[1:(k-1)])+j]*G1$S[[sum(n.S[1:(k-1)])+j]]
        }
      }
    }
    first <- G1$smooth[[k]]$first.para
    last <- G1$smooth[[k]]$last.para
    Lambda1[first:last, first:last] <- Lam[[k]]
  }
  
  n.smooth <- length(G2$smooth)
  Lambda2 <- matrix(0, np2, np2)
  Lam <- list(); n.S <- numeric(n.smooth) # penalty matrix
  for(k in 1:n.smooth) {
    n.S[k] <- length(G2$smooth[[k]]$S)
    if(k==1) {
      Lam[[k]] <- sp2[k]*G2$S[[k]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp2[j]*G2$S[[j]]
        }
      }
    }
    else {
      Lam[[k]] <- sp2[sum(n.S[1:(k-1)])+1]*G2$S[[sum(n.S[1:(k-1)])+1]]
      if(n.S[k]>1) {
        for(j in 2:n.S[k]) {
          Lam[[k]] <- Lam[[k]]+sp2[sum(n.S[1:(k-1)])+j]*G2$S[[sum(n.S[1:(k-1)])+j]]
        }
      }
    }
    first <- G2$smooth[[k]]$first.para
    last <- G2$smooth[[k]]$last.para
    Lambda2[first:last, first:last] <- Lam[[k]]
  }
  
  DS1 <- diag(eigen(Lambda1)$values[abs(eigen(Lambda1)$values)>1e-10])
  DS2 <- diag(eigen(Lambda2)$values[abs(eigen(Lambda2)$values)>1e-10])
  
  X1 <- G1$X; X2 <- G2$X
  
  loglik <- loglikfun(y, mu, p) # log-likelihood
  ploglik <- loglik - as.numeric(0.5*t(b1)%*%Lambda1%*%b1) -  as.numeric(0.5*t(b2)%*%Lambda2%*%b2)
  
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
  I.theta[(np1+1):(np1+np2),(np1+1):(np1+np2)] <- t(X2) %*% G.rho.p %*% X2 - Lambda2
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

