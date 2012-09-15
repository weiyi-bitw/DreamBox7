fwMCMC = function(x, surv, kmax = 2, delta = 0.01, maxIter = 10000){
	w = coxph(surv~., data=data.frame(t(x)))$coef
	w = w/w[1]
	m = nrow(x)
	n = ncol(x)
	out = .C("findWeightMCMC", x=as.double(x), surv=as.double(surv),w = as.double(w), mIn = as.integer(m), nIn = as.integer(n),kmaxIn=as.integer(kmax), dIn = as.double(delta), maxIterIn = as.integer(maxIter))
	return (out)
}
