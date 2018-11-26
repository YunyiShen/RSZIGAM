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

Hessian_sum_helper_pois = function(detmat,lambda,p,N){
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
		for(j in 1:ncol(detmat)){
			p.sumsqr[i,j] = (sum(fns * (detmat[i,j]-Ns*p[i,j])/(p[i,j]*(1-p[i,j]))))^2
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


occu.post.weight_helper = function(det.vec,p.vec,psi){
	log_detli = det.vec * log(p.vec) + (1-det.vec) * log(1-det.vec)
	detli = exp(log_detli)
	Zr = psi * detli + (1-psi) * (max(detvec)==0)
	return(psi * detli/Zr)
}
