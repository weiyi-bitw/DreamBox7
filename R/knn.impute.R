cwknn.impute <- function(X, k = 11){
	idx = which(is.na(X))
	nmiss = length(idx)
	idxx = idx %% nrow(X)
	idxx[idxx==0] = nrow(X)
	idxy = ceiling(idx / nrow(X))
	
	m = nrow(X)
	out = rep(NA, nmiss)	

	for(i in 1:nmiss){
		# cat(i, "/", nmiss, "\n");
		flush.console()
		ii = which(idxx == idxx[i])
		y_garbage = idxy[ii]
		db = X[, -y_garbage]
		vec = X[,idxy[i]]
		n = ncol(db)
		out[i] = .Call("cwknnR2C", db, vec, k, m, n, idxx[i]-1)
	}

	X[idx] = out
}

ewknn.impute <- function(X, wvec = rep(1, nrow(X)), k=11){
	idx = which(is.na(X))
	nmiss = length(idx)
	idxx = idx %% nrow(X)
	idxx[idxx==0] = nrow(X)
	idxy = ceiling(idx / nrow(X))
	
	m = nrow(X)
	out = rep(NA, nmiss)	

	for(i in 1:nmiss){
		# cat(i, "/", nmiss, "\n");
		flush.console()
		ii = which(idxx == idxx[i])
		y_garbage = idxy[ii]
		db = X[, -y_garbage]
		vec = X[,idxy[i]]
		n = ncol(db)
		out[i] = .Call("ewknnR2C", db, vec, wvec, k, m, n, idxx[i]-1)
	}

	X[idx] = out
}

ewknn.predict <- function(X, t, qX, wvec=rep(1, nrow(X)), k=11){
	
	m = nrow(X)
	n = ncol(X)
	
	if(length(t) != n){stop("error: length(t) != ncol(X)")}
	if(nrow(qX) != m){stop("error: nrow(qX) != nrow(X)!!")}
	if(length(wvec) != m){stop("error: length(wvec) != nrow(X)")}
	
	nq = ncol(qX)
	out = .Call("ewknnPredictR2C", X, t, qX, wvec, k, m, n, nq)
	
	return(out)
}

preproClncKNN=function(c, survobj){
	n = ncol(c)
	m = nrow(c)
	isFactor = rep(NA, n)
	for(i in 1:n) isFactor[i] = is.factor(c[,i])
	names(isFactor) = colnames(c)
	isFactor["grade"] = TRUE
	c[,"grade"] = factor(c[,"grade"])
	distWeight = list()
	ws = rep(NA, n)
	names(ws) = colnames(c)
	for(i in 1:n){
		x = c[,i]
		if(isFactor[i]){
			max.i = which.max(table(x))
			levx = levels(x)
			contrmat = matrix(0, nrow=length(levx), ncol=length(levx)-1)
			contrmat[max.i,] = rep(-1, ncol(contrmat))
			contrmat[-max.i,] = diag(1, ncol(contrmat-1))
			contrasts(x) = contrmat
			fit = coxph(survobj~x)
			coeff = fit$coef
			coeff[is.na(coeff)] = 0
			coeffFac = rep(0, length(levx))
			coeffFac[-max.i] = coeff
			coeffFac[max.i] = 0-sum(coeff)
			coeff = coeffFac
			names(coeff) = levels(x) 
		}else{
			fit = coxph(survobj~x)
			coeff = fit$coef
			coeff[is.na(coeff)] = 0
		}
		ws[i] = survConcordance(survobj~predict(fit))$concordance
		distWeight[[i]] = coeff
	}
	names(distWeight) = colnames(c)
	out = list(clinical = c, isFactor = isFactor, distWeight = distWeight, weights = ws)
}

clncDist = function(a, b, isFac, dw){
	
	d = mapply(function(x, y, isFaci, dwi){
		if(isFaci) abs(dwi[x] - dwi[y])
		else abs(x - y) * dwi
	},
	a, b, isFac, dw
	)
	mean(d^2)
}
