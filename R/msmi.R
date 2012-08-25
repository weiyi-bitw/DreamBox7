getMSMI = function(x, survobj){
	#t = survobj[,1]
	#s = survobj[,2]
	n = length(x)
	if (nrow(survobj) != n){
		stop("Two vectors have different length!")
	}
	tstart = sort(survobj[,1])[20]
	tend = sort(survobj[, 1], decreasing=T)[20]
	ti = tstart	
	mi = 0
	tmax = 0
	while(ti < tend){
		idx = which(survobj[,1] <= ti & survobj[,2]==0)
		b = rep(1, n)
		b[survobj[,1]<= ti] = 0
		mii = getMI(x[-idx], b[-idx], normalize=T)
		if (mii > mi){
			mi = mii
			tmax = ti
		}
		ti = ti + 0.1
	}
	out = list(mi = mi, t.max = tmax)
	return (out)

}
