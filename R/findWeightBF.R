#=========================================
#
#  BFFW - Brute-Force Finding Weights
#
#=========================================

BFFW = function(x, surv, w=NULL, delta = 0.5, maxIter = 10000, verbose = TRUE){
	if(is.null(w)){
		w = coxph(surv~., data=data.frame(t(x)))$coef
		w = w/w[1]
	}
	m = nrow(x)
	n = ncol(x)
	out = .C("findWeightBF", x=as.double(x), surv=as.double(surv),w = as.double(w), mIn = as.integer(m), nIn = as.integer(n), dIn = as.double(delta), maxIterIn = as.integer(maxIter), verbose = as.integer(verbose))
	return (out$w)
}

#=========================================
#
# BFFS - Brute-Force Feature Selection
#
#=========================================

BFFS = function(x, surv, folds = 5){
	n = nrow(x)
	m = ncol(x)
	
	idx = sample(1:n, floor(4*n/5))
	ccdi = getAllCCDIWz(t(x[idx,]), surv[idx,])
	
	bestFeatures = which.max(ccdi)
	bestNumFeatures = 1
	bestCCDI = max(ccdi)
	
	numFeatures = 2
	

	while(TRUE){
		currentBestCCDI = 0
		currentBestAddFeature = 0
		for(i in 1:m){
			if(! (i %in% bestFeatures) ){
				ft = c(bestFeatures, i)
				cm = coxph(surv[idx,]~., data=x[idx, ft])
				pred = predict(cm, x[-idx,ft])
				ccdi = getCCDIdx(pred, surv[-idx,])
				if(ccdi > currentBestCCDI){
					currentBestCCDI = ccdi
					currentBestAddFeature = i
				}
			}
		}
		if(currentBestCCDI > bestCCDI){
			bestFeatures = c(bestFeatures, currentBestAddFeature)
			bestNumFeatures = bestNumFeatures + 1
			bestCCDI = currentBestCCDI
			cat("New Best Feature:\n")
			cat(colnames(x)[bestFeatures], "\n")
			cat("New Best CCDI:\n")
			cat(bestCCDI, "\n");flush.console()
		}else{
			break
		}


	}
	cat("Done.\n");flush.console()
	list(bestFeatures, bestNumFeatures, bestCCDI)
}
