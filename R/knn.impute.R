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
