cwknn.impute <- function(X, k = 10){
	idx = which(is.na(X))
	nmiss = length(idx)
	idxx = idx %% nrow(X)
	idxx[idxx==0] = nrow(X)
	idxy = ceiling(idx / nrow(X))
	
	m = nrow(X)
	out = rep(NA, nmiss)	

	for(i in 1:nmiss){
		cat(i, "/", nmiss, "\n");
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
