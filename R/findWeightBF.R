#=========================================
#
#  BFFW - Brute-Force Finding Weights
#
#=========================================

BFFW = function(x, surv, w=NULL, delta = 0.5, maxIter = 10000, verbose = TRUE){
	if(is.null(w)){
		w = coxph(surv~., data=data.frame(t(x)))$coef
		w[is.na(w)] = 1
	}
	mid = order(abs(w))[ceiling(length(w)/2)]
	w = w/abs(w[mid])
	temp = x[mid,]
	x[mid,] = x[1,]
	x[1,] = temp

	temp = w[mid]
	w[mid] = w[1]
	w[1] = temp
	
	m = nrow(x)
	n = ncol(x)
	out = .C("findWeightBF", x=as.double(x), surv=as.double(surv),w = as.double(w), mIn = as.integer(m), nIn = as.integer(n), dIn = as.double(delta), maxIterIn = as.integer(maxIter), verbose = as.integer(verbose))

	wo = out$w
	temp = wo[1]
	wo[1] = wo[mid]
	wo[mid] = temp

	return (wo)
}

BFFWFS = function(x, surv, numFeatures = 5, randomShuffle = 10000, maxIter=1000, k = 2, verbose=TRUE, ft.init= NULL,seed=913){
	set.seed(seed)
	m = nrow(x) # number of samples
	n = ncol(x) # number of features

	out = list()

	#idx.train = sample(1:n, floor(train.fraction * n))
	count = 0
	if(is.null(ft.init)){
		ft = 1:numFeatures
	}else{
		ft = ft.init
	}
	#cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
	#p = predict(cm, x[-idx.train,ft])
	#ccdi = getCCDIdx(p, surv[-idx.train,])
	w = BFFW(x[ft,], surv, maxIter = maxIter, verbose=FALSE)
	ccdi = getCCDIdx(w %*% x[ft,], surv)
	cat("Starting CCDI:", ccdi, "\n");flush.console()
	bestCCDI = ccdi
	bestFeature = ft
	bestW = w

	while( count < randomShuffle ){
		idx = sample(1:numFeatures, k)
		ft[idx] = sample(setdiff(1:m, ft[-idx]), k)
		#cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
		#p = predict(cm, x[-idx.train,ft])
		#ccdi = getCCDIdx(p, surv[-idx.train,])
		w = BFFW(x[ft,], surv, maxIter = maxIter, verbose=FALSE)
		ccdi = getCCDIdx(w %*% x[ft,], surv)
		if(ccdi > bestCCDI){
			cat("New CCDI:\t", ccdi, "\n");flush.console()
			bestCCDI = ccdi
			cat("New feature: \n\t", rownames(x)[ft], "\n");flush.console()
			bestFeature = ft
			bestCM = w
			count = 0
		}else{
			ft = bestFeature
		}
		count = count + 1
		if(count %% 1000 == 0) cat(count, "\n");flush.console()
	}

	out[[r]] = list(bestFeatures=rownames(x)[bestFeature], bestCCDI=bestCCDI, bestCM = bestCM)
	cat("Done.\n");flush.console()
	return (out)
}
#=========================================
#
# BFFS - Brute-Force Feature Selection
#
#=========================================

BFFS = function(x, surv, numFeatures = 5, randomShuffle = 10000, k = 2, verbose=TRUE, train.fraction=0.5, rounds = 10, idx.train = NULL, ft.init=NULL,  seed=913){
	set.seed(seed)
	n = nrow(x) # number of samples
	m = ncol(x) # number of features

	out = list()

	if(is.null(idx.train)){


	for(r in 1:rounds){
	#idx.train = sample(1:n, floor(train.fraction * n))
	count = 0
	if(is.null(ft.init)){
		ft = 1:numFeatures
	}else{
		ft = ft.init
	}
	#cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
	#p = predict(cm, x[-idx.train,ft])
	#ccdi = getCCDIdx(p, surv[-idx.train,])
	cm = coxph(surv~., data=x[,ft])
	ccd = cm$concordance
	ccdi = (ccd[1]) / (ccd[1] + ccd[2])
	bestCCDI = ccdi
	bestFeature = ft
	bestCM = cm

	while( count < randomShuffle ){
		idx = sample(1:numFeatures, k)
		ft[idx] = sample(setdiff(1:m, ft[-idx]), k)
		#cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
		#p = predict(cm, x[-idx.train,ft])
		#ccdi = getCCDIdx(p, surv[-idx.train,])
		cm = coxph(surv~., data=x[,ft])
		ccd = cm$concordance
		ccdi = ccd[1] / (ccd[1] + ccd[2])
		if(ccdi > bestCCDI){
			cat("New CCDI:\t", ccdi, "\n");flush.console()
			bestCCDI = ccdi
			cat("New feature: \n\t", colnames(x)[ft], "\n");flush.console()
			bestFeature = ft
			bestCM = cm
			count = 0
		}else{
			ft = bestFeature
		}
		count = count + 1
		if(count %% 1000 == 0) cat(count, "\n");flush.console()
	}

	out[[r]] = list(bestFeatures=colnames(x)[bestFeature], bestCCDI=bestCCDI, bestCM = bestCM)

	}
	cat("Done.\n");flush.console()
	return (out)


	}else{
		count = 0
		ft = 1:numFeatures
		cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
		p = predict(cm, x[-idx.train,ft])
		ccdi = getCCDIdx(p, surv[-idx.train,])
		bestCCDI = ccdi
		bestFeature = ft
		bestCM = cm

		while( count < randomShuffle ){
			idx = sample(1:numFeatures, k)
			ft[idx] = sample(setdiff(1:m, ft[-idx]), k)
			cm = coxph(surv[idx.train,]~., data=x[idx.train,ft])
			p = predict(cm, x[-idx.train,ft])
			ccdi = getCCDIdx(p, surv[-idx.train,])
			if(ccdi > bestCCDI){
				cat("New CCDI:\t", ccdi, "\n");flush.console()
				bestCCDI = ccdi
				cat("New feature: \n\t", colnames(x)[ft], "\n");flush.console()
				bestFeature = ft
				bestCM = cm
				count = 0
			}else{
				ft = bestFeature
			}
		count = count + 1
		if(count %% 1000 == 0) cat(count, "\n");flush.console()
		}

		out = list(bestFeatures=colnames(x)[bestFeature], bestCCDI=bestCCDI, bestCM = bestCM)

		cat("Done.\n");flush.console()
		return (out)

	}
}

BICFS = function(x, surv, verbose=FALSE, rounds = 10, train.fraction=0.4, seed=121229, k = 2){
	n = nrow(x)
	m = ncol(x)

	foldIndices = list()
	for(i in 1:rounds){
		set.seed(seed + i)
		foldIndices[[i]] = sample(1:n, floor(n*train.fraction))
	}
	out = list()
	for(r in 1:rounds){
		if(verbose) cat("Round" , r , '...\n');flush.console()
		idx.train=foldIndices[[r]]
		if(verbose) print(idx.train)
		xt = x[idx.train,]
		survt = surv[idx.train,]
		upper = terms(survt~., data=xt)
		cm = step(coxph(survt~1, data=xt), scope=upper, direction="both", k=k, trace=verbose)
		out[[r]] = list(bestFeatures = attr(cm$term, "term.labels"), bestCM = cm)
	}
	if(verbose)cat("Done.\n");flush.console()
	return (out)
}


#=========================================
#
# qBFFS - Brute-Force Feature Selection cluster version
#
#=========================================

qBFFS = function(x, surv, folds = 10, randomShuffle = 10000, kmax = 2, verbose=TRUE, thisSeg=0, totalSeg=1, seed=92){
	n = nrow(x) # number of samples
	m = ncol(x) # number of features
	
    	set.seed(seed)
	foldIndices = createFolds(1:n, k=folds, list=TRUE)
	
	start = floor(folds * thisSeg / totalSeg)
	end = floor(folds * (thisSeg+1) / totalSeg)
	for(f in 1:folds){
		if(f <= start | f > end) next
		idx = foldIndices[[f]]
		cat("======= Fold:", f, "/", folds, "=======\n");flush.console()
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
		filename = paste("bestFeatures.",  sprintf("%02d", as.numeric(f) ), ".txt", sep="" )
		write(bestFeatures, file=filename, sep='\t')
		write(bestCCDI, file=filename)
		write(bestTestCCDI, file=filename)
		unlink(filename)

		cat("Done.\n");flush.console()
		

	}
}
