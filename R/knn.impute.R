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

ewknn.predict <- function(X, t, qX, wvec=rep(1, nrow(X)), k=ncol(X)){
	
	m = nrow(X)
	n = ncol(X)
	
	if(length(t) != n){stop("error: length(t) != ncol(X)")}
	if(nrow(qX) != m){stop("error: nrow(qX) != nrow(X)!!")}
	if(length(wvec) != m){stop("error: length(wvec) != nrow(X)")}
	
	nq = ncol(qX)
	out = .Call("ewknnPredictR2C", X, t, qX, wvec, k, m, n, nq)
	
	return(out)
}


preproClncKNN=function(c, survobj=NULL, isFactorIn=NULL, dwIn=NULL){
	n = ncol(c)
	m = nrow(c)
	if(!is.null(survobj)) if(nrow(survobj) != m) {stop("error: nrow(survobj) != nrow(clnc)")}
	if(!is.null(dwIn)) if(length(dwIn) != n) stop("error: length(distWeight) != ncol(clnc)")
	if(!is.null(isFactorIn)) if(length(isFactorIn) != n) stop("error: length(isFactor) != ncol(clnc)")
	
	if(is.null(survobj) & is.null(dwIn)) stop("error: distance weight and survival object cannot be both NULL!")
	if(is.null(isFactorIn)){
		isFactor = rep(NA, n)
		for(i in 1:n) isFactor[i] = is.factor(c[,i])
		names(isFactor) = colnames(c)
	}else{
		isFactor = isFactorIn
	}
	distWeight = list()
	ccdi = rep(NA, n)
	names(ccdi) = colnames(c)
	for(i in 1:n){
		#cat(colnames(c)[i], "\n");flush.console()
		x = c[,i]
		if(isFactor[i]){
			max.i = which.max(table(x))
			levx = levels(x)
			contrmat = matrix(0, nrow=length(levx), ncol=length(levx)-1)
			contrmat[max.i,] = rep(-1, ncol(contrmat))
			contrmat[-max.i,] = diag(1, ncol(contrmat-1))
			contrasts(x) = contrmat
			if(!is.null(survobj)){
				fit = coxph(survobj~x)
				coeff = fit$coef
				coeff[is.na(coeff)] = 0
				coeffFac = rep(0, length(levx))
				coeffFac[-max.i] = coeff
				coeffFac[max.i] = 0-sum(coeff)
				coeff = coeffFac
				names(coeff) = levels(x)
			}else{
				coeff = dwIn[[colnames(c)[i]]]
			}
			c[,i] = coeff[c[,i]]
		}else{
			if(!is.null(survobj)){
				fit = coxph(survobj~x)
				coeff = fit$coef
			}else{
				coeff = dwIn[[colnames(c)[i]]]
			}
			c[,i] = c[,i] * coeff
		}
		if(!is.null(survobj)) ccdi[i] = survConcordance(survobj~predict(fit))$concordance
		distWeight[[i]] = coeff
	}
	names(distWeight) = colnames(c)
	if(!is.null(dwIn)) distWeight = dwIn
	if(!is.null(survobj)) {
		idx = which(survobj[,2]==1)
		cout = c[idx,]
		tout = survobj[idx, 1]
		idx2 = which(ccdi > 0.55 | ccdi < 0.45)
		cout = cout[, names(idx2)]
		ccdi = ccdi[names(idx2)]
		isFactor = isFactor[names(idx2)]
		distWeight = distWeight[names(idx2)]
	}else{
		cout = c
		tout = NULL
	}
	
	out = list(clinical = cout, time = tout, concordance=ccdi, isFactor = isFactor, distWeight = distWeight)
}

clncDist = function(a, b, isFac, dw, wvec, minimumDist=0.0001){
	
	d = mapply(function(x, y, isFaci, dwi, wvec){
		if(isFaci) abs(dwi[x] - dwi[y])*wvec
		else abs(as.numeric(x) - as.numeric(y)) * dwi *wvec
	},
	a, b, isFac, dw, wvec
	)
	sum(wvec)/sqrt(sum(d^2 + minimumDist))
}

getAllClncDist = function(C, qc, isFactor, distWeights, wvec){
	apply( C, 1, function(c, qc, isFac, dw, wv){
		clncDist(c, qc, isFac, dw, wv) 
	}, qc = qc, isFac = isFactor, dw = distWeights, wv = wvec
	)
}

ewknn.predict.clnc <- function(preproClnc, qC, wvec=rep(1, length(preproClnc$distWeight)), k=nrow(preproClnc$clinical)){
	clnc = t(preproClnc$clinical)
	m = nrow(clnc)
	n = ncol(clnc)
	
	qC = qC[, rownames(clnc)]

	if(ncol(qC) != m){stop("error: ncol(qC) != ncol(clnc)!!")}
	
	qC = t(preproClncKNN(qC, isFactorIn = preproClnc$isFactor, dwIn = preproClnc$distWeight)$clinical)
	nq = ncol(qC)
	t = .Call("ewknnPredictR2C", clnc, preproClnc$time, qC,wvec, k, m, n, nq)	

	return (t)	
}

