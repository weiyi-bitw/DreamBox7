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

BFFS = function(x, surv, folds = 10, randomShuffle = 10000, kmax = 2, verbose=TRUE){
	n = nrow(x) # number of samples
	m = ncol(x) # number of features
	
	foldIndices = createFolds(1:n, k=folds, list=TRUE)

	CRResults = foreach(idx = foldIndices) %do% {

		bestFeatures = 1
		bestNumFeatures = 1
		bestCCDI = 0
		bestTestCCDI = 0
		bestW = NULL
		while(TRUE){
			cat("Test", bestNumFeatures+1, "features...\n");flush.console()
			currentBestCCDI = 0
			currentBestFeature = NULL
			currentBestW = NULL
			currentBestTestCCDI = 0
			currentBestModel = NULL
			if(choose(m, bestNumFeatures+1) < randomShuffle){
				allComb = combn(1:m, bestNumFeatures+1)
				allCCDI = apply(allComb, 2, function(ii){
					cm = coxph(surv[-idx,]~., data=x[-idx, ii])
					#pred = predict(cm, x[idx, ii])
					#w = cm$coeff
					#w = BFFW(x[ii, -idx], surv[-idx,], verbose = FALSE)
					#pred = w %*% x[ii, idx]
					#getCCDIdx(pred, surv[idx,])
					ccd = cm$concordance
					(ccd[1] + ccd[3]/2) / (ccd[1] + ccd[2] + ccd[3])
				})
				best = which.max(allCCDI)
				currentBestCCDI = max(allCCDI)
				currentBestFeature = allComb[,best]
				currentBestModel = coxph(surv[-idx,]~., data=x[-idx, currentBestFeature])
				pred = predict(currentBestModel, x[idx,currentBestFeature])
				currentBestTestCCDI = getCCDIdx(pred, surv[idx,])
				currentBestW = currentBestModel$coeff
			}else{
				cat("Finding new features incrementally ... ");flush.console()
				for(i in 1:m){
					if(! (i %in% bestFeatures) ){
						ft = c(bestFeatures, i)
						cm = coxph(surv[-idx,]~., data=x[-idx, ft])
						#pred = predict(cm, x[idx,ft])
						w = cm$coeff
						#w = BFFW(x[ft,-idx], surv[-idx,], verbose=FALSE)
						#pred = w %*% x[ft,idx]
						#ccdi = getCCDIdx(pred, surv[idx,])
						ccd = cm$concordance
						ccdi = (ccd[1] + ccd[3]/2) / (ccd[1] + ccd[2] + ccd[3])
						if(ccdi > bestCCDI){
							pred = predict(cm, x[idx,ft])
							testccdi = getCCDIdx(pred, surv[idx,])
							if(testccdi > currentBestTestCCDI){
								currentBestCCDI = ccdi
								currentBestFeature = ft
								currentBestW = w
								currentBestModel = cm
								currentBestTestCCDI = testccdi
							}
						}
					}	
				}
				cat("Randomly shuffle...\n");flush.console()
				k = min(kmax, ceiling(bestNumFeatures/2))
				for(i in 1:randomShuffle){
					ft = currentBestFeature
					ft[sample(1:(bestNumFeatures+1), k)] = sample(setdiff(1:m, currentBestFeature), k)
					cm = coxph(surv[-idx,]~., data=x[-idx, ft])
					#pred = predict(cm, x[idx,ft])
					w = cm$coeff
					#w = BFFW(x[ft, -idx], surv[-idx,], verbose=FALSE)
					#pred = w %*% x[ft, idx]
					#ccdi = getCCDIdx(pred, surv[idx,])
					ccd = cm$concordance
					ccdi = (ccd[1] + ccd[3]/2) / (ccd[1] + ccd[2] + ccd[3])
					if(verbose) cat("Features:", colnames(x)[ft], "\tCCDI:", ccdi, "\n");flush.console()
					
					if(ccdi > bestCCDI){
						pred = predict(cm, x[idx,ft])
						testccdi = getCCDIdx(pred, surv[idx,])
						if(verbose) cat("\tTest CCDI:", testccdi, "\n");flush.console()
						if(testccdi > currentBestTestCCDI){
							currentBestCCDI = ccdi
							currentBestFeature = ft
							currentBestW = w
							currentBestModel = cm
							currentBestTestCCDI = testccdi
						}
					}
				}
			}
			#cat(currentBestCCDI, '\t', bestCCDI, '\n');flush.console()
			cat("New Best Feature:\n")
			cat(colnames(x)[currentBestFeature], "\n")
			cat("New Best CCDI:\n")
			cat(currentBestCCDI, " / ", currentBestTestCCDI, "\n");flush.console()
			if(currentBestTestCCDI > bestTestCCDI){
				bestFeatures = currentBestFeature
				bestNumFeatures = bestNumFeatures + 1
				bestCCDI = currentBestCCDI
				bestTestCCDI = currentBestTestCCDI
				bestW = currentBestW
			}else{
				break
			}

		}	
		cat("Done.\n");flush.console()
		return (list(bestFeatures=colnames(x)[bestFeatures], bestNumFeatures=bestNumFeatures, bestCCDI=bestCCDI, bestTestCCDI = bestTestCCDI, bestW = bestW))
	}

	CRFeatures = sapply( CRResults, function(x) x$bestFeatures )
	CRNumFeatures = sapply( CRResults, function(x) x$bestNumFeatures )
	CRCCDI = sapply(CRResults, function(x) x$bestCCDI )
	CRTestCCDI = sapply(CRResults, function(x) x$bestTestCCDI)
	CRBestW = sapply(CRResults, function(x) x$bestW )

	return(list(CRFeatures = CRFeatures, CRNumFeatures=CRNumFeatures, CRCCDI = CRCCDI, CRTestCCDI = CRTestCCDI, CRBestW = CRBestW, foldIndices = foldIndices))
}
