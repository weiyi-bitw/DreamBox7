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

BFFS = function(x, surv, folds = 10, randomShuffle = 10000){
	n = nrow(x) # number of samples
	m = ncol(x) # number of features
	
	foldIndices = createFolds(1:n, k=folds, list=TRUE)

	CRResults = foreach(idx = foldIndices) %do% {

		bestFeatures = 1
		bestNumFeatures = 1
		bestCCDI = 0
		bestW = NULL
		while(TRUE){
			cat("Test", bestNumFeatures+1, "features...\n");flush.console()
			currentBestCCDI = 0
			currentBestFeature = NULL
			currentBestW = NULL
			if(choose(m, bestNumFeatures+1) < randomShuffle){
				allComb = combn(1:m, bestNumFeatures+1)
				allCCDI = apply(allComb, 2, function(ii){
					cm = coxph(surv[-idx,]~., data=x[-idx, ii])
					pred = predict(cm, x[idx, ii])
					w = cm$coeff
					#w = BFFW(x[ii, -idx], surv[-idx,], verbose = FALSE)
					#pred = w %*% x[ii, idx]
					getCCDIdx(pred, surv[idx,])
				})
				best = which.max(allCCDI)
				currentBestCCDI = max(allCCDI)
				currentBestFeature = allComb[,best]
			}else{
				cat("Finding new features incrementally ... ");flush.console()
				for(i in 1:m){
					if(! (i %in% bestFeatures) ){
						ft = c(bestFeatures, i)
						cm = coxph(surv[-idx,]~., data=x[-idx, ft])
						pred = predict(cm, x[idx,ft])
						w = cm$coeff
						#w = BFFW(x[ft,-idx], surv[-idx,], verbose=FALSE)
						#pred = w %*% x[ft,idx]
						ccdi = getCCDIdx(pred, surv[idx,])
						if(ccdi > currentBestCCDI){
							currentBestCCDI = ccdi
							currentBestFeature = ft
							currentBestW = w
						}
					}	
				}
				cat("Randomly shuffle...\n");flush.console()
				k = min(3, ceiling(bestNumFeatures/2))
				for(i in 1:randomShuffle){
					ft = currentBestFeature
					ft[sample(1:(bestNumFeatures+1), k)] = sample(setdiff(1:m, currentBestFeature), k)
					cm = coxph(surv[-idx,]~., data=x[-idx, ft])
					pred = predict(cm, x[idx,ft])
					w = cm$coeff
					#w = BFFW(x[ft, -idx], surv[-idx,], verbose=FALSE)
					#pred = w %*% x[ft, idx]
					ccdi = getCCDIdx(pred, surv[idx,])
					if(ccdi > currentBestCCDI){
						currentBestCCDI = ccdi
						currentBestFeature = ft
						currentBestW = w
					}
				}
			}
			#cat(currentBestCCDI, '\t', bestCCDI, '\n');flush.console()
			if(currentBestCCDI > bestCCDI){
				bestFeatures = currentBestFeature
				bestNumFeatures = bestNumFeatures + 1
				bestCCDI = currentBestCCDI
				bestW = currentBestW
				cat("New Best Feature:\n")
				cat(colnames(x)[bestFeatures], "\n")
				cat("New Best CCDI:\n")
				cat(bestCCDI, "\n");flush.console()
			}else{
				break
			}

		}	
		cat("Done.\n");flush.console()
		return (list(bestFeatures=colnames(x)[bestFeatures], bestNumFeatures=bestNumFeatures, bestCCDI=bestCCDI, bestW = bestW))
	}

	CRFeatures = sapply( CRResults, function(x) x$bestFeatures )
	CRNumFeatures = sapply( CRResults, function(x) x$bestNumFeatures )
	CRCCDI = sapply(CRResults, function(x) x$bestCCDI )
	CRBestW = sapply(CRResults, function(x) x$bestW )

	return(list(CRFeatures = CRFeatures, CRNumFeatures=CRNumFeatures, CRCCDI = CRCCDI, CRBestW = CRBestW))
}
