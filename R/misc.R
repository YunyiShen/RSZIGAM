likelihood.fnr = function(n,det.vec,lambda,p.vec){
	det.p = prod(dbinom(det.vec,rep(n,length(det.vec)),p.vec))
	pois  = dpois(n,lambda)
	return(det.p * pois)
}

post.weight_helper = function(n,det.vec, lambda,p.vec,psi,N){
	Ns = as.matrix( max(det.vec):N )
	fns = apply(Ns,1,likelihood.fnr,det.vec,lambda,p.vec)
	Zr = psi*sum(fns) + (1-psi) * (sum(det.vec)==0)
	fn = likelihood.fnr(n,det.vec,lambda,p.vec)
	return(psi*fn/Zr)
}

log_likelihood_pois_eachsite = function(det.vec,lambda,p.vec,psi,N){
	# n.site = nrow(detmat)
	nvec = max(det.vec):N
	gr = apply(as.matrix(nvec),1,likelihood.fnr,det.vec=det.vec,lambda=lambda,p.vec = p.vec)
	gr = sum(gr) # given all N the probability of having data, which is actually N-mixture
	e = 1.0*(max(det.vec)!=0)
	logL = e * (log(psi) + log(gr))+ (1-e)*(log(1-psi + psi * gr)) # log likelihood with zero inflating 
	return(logL)
}

log_likelihood_pois = function(detmat,lambda,p,psi,N){
	n.site = nrow(detmat)
	logL = matrix(0,n.site,1)
	for(i in 1:n.site){
		logL[i] = log_likelihood_pois_eachsite(det.vec = detmat[i,],lambda=lambda,p.vec = p[i,],psi = psi[i],N=N) 
	}
	return(logL)
}

Hessian_sum_helper_pois = function(detmat,lambda,p,N){# Useful sum in calculating Hessians, before really calculating it. Actually we can know that design matrix X involved. And Hessian has qauadratic form of X.
	n.site = nrow(detmat)
	lambda.sqrsumfnn = 0*detmat[,1]
	lambda.sumsqrfnn = lambda.sqrsumfnn
	psi.fnminusId0 = lambda.sqrsumfnn
	p.sumsqr=0 * detmat
	p.sum = 0 * detmat
	
	for(i in 1:n.site){ # this useful matrix has n.site numbers for psi and lambda, n.site*n.period for p
		Ns = as.matrix(max(detmat[i,]):N)
		fns = apply(Ns,1,likelihood.fnr,det.vec=detmat[i,],lambda[i],p.vec=p[i,])
		nminusmu = Ns-lambda[i]
		lambda.sqrsumfnn[i] = (sum(fns*nminusmu))^2
		lambda.sumsqrfnn[i] = sum(fns * nminusmu^2)
		psi.fnminusId0[i] = sum(fns)-(max(detmat[i,])==0)
		for(j in 1:ncol(detmat)){
			p.sumsqr[i,j] = (sum(fns * (detmat[i,j]-Ns*p[i,j])/(p[i,j]*(1-p[i,j]))))^2 # see document for derivation
			p.sum[i,j] = (1/(p[i,j]*(1-p[i,j]))^2) * sum(fns * ((detmat[i,j]-Ns*p[i,j])^2-Ns*(1-p[i,j]*p[i,j]-(1-2*p[i,j])*(detmat[i,j]-Ns * p[i,j]))))
		}
	}
	res = list(lambda.sqrsumfnn,lambda.sumsqrfnn,psi.fnminusId0,p.sumsqr,p.sum)
	names(res) = c('lambda.sqrsumfnn','lambda.sumsqrfnn','psi.fnminusId0','p.sumsqr','p.sum')
	return(res)
}

# helper functions for occupancy
log_likelihood_occu = function(detmat,p,psi){
	log_detli = det.vec * log(p.vec) + (1-det.vec) * log(1-det.vec)
	detli = exp(log_detli)
	Zr = psi * detli + (1-psi) * (max(detvec)==0)
	return(log(Zr))
}

log_likelihood_occu_eachsite = function(det.vec,p.vec,psi){
	# n.site = nrow(detmat)
	log_gr =exp( det.vec*log(p.vec) + (1-det.vec)*log(1-p.vec))
	e = 1.0*(max(det.vec)!=0)
	logL = e * (log(psi) + log(gr))+ (1-e)*(log(1-psi + psi * gr)) # log likelihood with zero inflating 
	return(logL)
}

Hessian_sum_helper_occu = function(detmat,p){
	allp = exp( rowSums( log(p)*detmat + log(1-p) * (1-detmat) ))

	res = list(allpminusId0=allp-(rowSums(detmat)==0),allp = allp)
	return(res)
}

occu.post.weight_helper = function(det.vec,p.vec,psi){
	log_detli = det.vec * log(p.vec) + (1-det.vec) * log(1-det.vec)
	detli = exp(log_detli)
	Zr = psi * detli + (1-psi) * (max(detvec)==0)
	return(psi * detli/Zr)
}

# other helpers 
meshgrid = function (xrange, yrange){ # meshgrid from matlab
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

check.data = function(data){
	cat("Checking data formation...\n")
	msg = "Data formation all pass.\n"
	allright = TRUE
	if(is.null(data$detmat)){
		allright = FALSE
		msg = "Missing detection matrix.\n"
		return(list(allright = allright,msg=msg))
	}
	
	if(is.null(data$envX)){
		allright = FALSE
		msg = "Missing environmental data.\n"
		return(list(allright=allright,msg=msg))
	}
	
	n.site = nrow(data$detmat)
	n.period = ncol(data$detmat)
	len.data = length(data)
	
	if(nrow(data$envX)!=n.site){
		allright = FALSE
		msg = "Environmental variables show site number which does not agree with detection matrix.\n"
		return(list(allright = allright,msg = msg))
	}
	
	if(len.data != (n.period + 2) & len.data != 2){
		allright = FALSE
		msg = "Detection variables show number of detection periods which does not agree with detection matrix.\n"
		return(list(allright = allright,msg = msg))
	}
	
	n.site.detvar = sapply(data[c(-1,-2)],function(mat,n.site){nrow(mat)!=n.site},n.site=n.site)
	if(sum(n.site.detvar)!=0){
		disagreeperiod = which(n.site.detvar)
		msg = paste("Detection variables of period",paste(disagreeperiod,collapse=","),"show site numbers which do not agree with detection matrix.\n")
		return(list(allright = allright,msg = msg))
	}
	
	return(list(allright = allright,msg = msg))
}