fwMCMC = function(x, surv, delta = 0.01, maxIter = 10000){
	w = coxph(surv~., data=data.frame(t(x)))$coef
	m = nrow(x)
	n = ncol(x)
	out = .C("findWeightMCMC", x=x, surv=surv,w = w, mIn = m, nIn = n, dIn = delta, maxIterIn = maxIter)
	return (out)
}
