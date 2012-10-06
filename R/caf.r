attractorScanningGL = function(data, genome, alpha=(2:12)/2, windowSize = 50, maxIter = 100, epsilon=1E-14, bin=6, so=3, score.position=5, num.output=10, negateMI=TRUE, verbose=TRUE, saveAllScore=FALSE){
	genes.genome = rownames(genome)
	data = data[intersect(genes.genome, rownames(data)),] # sort the rownames according to genomic location
	m = nrow(data)
	n = ncol(data)

	genes.data = rownames(data)
	mg = nrow(genome)
	alist = list()
	ascore = rep(0, m)
	aalpha = rep(NA, m)
	for(i in 1:m){
		if(verbose & i %% 100 == 0){
			cat(i , "/", m, "\n"); flush.console()
		}
		idx.center = which(genes.genome == genes.data[i])
		winspan = floor(windowSize/2)
		wd = genes.genome[max(1, idx.center - winspan) : min(mg, idx.center+winspan)]
		wd = intersect(wd, genes.data)
		if(length(wd) < num.output){
			alist[[i]] = NA
			next
		}
		dataIn = data[wd,]
		mi = getAllMIWz(dataIn, dataIn[genes.data[i],], bin=bin, so=so, negateMI=negateMI)
		premi = mi
		for(a in alpha){
			w = abs(mi)^a / sum(abs(mi)^a)
			w[mi < 0] = 0
			metagene = w %*% dataIn
			c = 0
			while(c < maxIter){
				mi = getAllMIWz(dataIn, metagene, bin=bin, so=so, negateMI=negateMI)
				delta = sum((mi - premi)^2)
				if(delta < epsilon){
					break
				}
				premi = mi
				mi[mi < 0] = 0
				metagene = w %*% dataIn
				c = c + 1
			}
			if (c >= maxIter){ 
				alist[[i]] = NA
			}else{
				aout = sort(mi, decreasing=T)[1:num.output]
				alist[[i]] = aout
				score = aout[score.position]
				if(score >= ascore[i]){
					ascore[i] = score
					aalpha[i] = a
				}
			}
		}
	}
	killidx = which(is.na(alist))
	if (length(killidx) > 0){
		alist[killidx] = NULL
		ascore = ascore[-killidx]
		aalpha = aalpha[-killidx]
		aglist = sapply(alist, names)
		colnames(aglist) = genes.data[-killidx]
		if (saveAllScore){
			 alist = simplify2array(alist)
			colnames(alist) = genes.data[-killidx]
		}else{
			alist=NULL
		}
	}else{
		aglist = sapply(alist, names)
		colnames(aglist) = genes.data
		if (saveAllScore){
			 alist = simplify2array(alist)
			colnames(alist) = genes.data
		}else{
			alist=NULL
		}

	}	
	
	return(list(attractome=aglist, score=ascore, bestAlphas = aalpha, scoremat = alist))

}


findAttractor = function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE){
	dataIn = data
	if(rankBased){
		vec = rank(vec)
		for(i in 1:nrow(data)){
			dataIn[i,] = rank(data[i,])
		} 
	}
	mi = getAllMIWz(dataIn, vec, bin=bin, so=so, negateMI = negateMI)
	premi = mi
	w = abs(mi)^a / sum(abs(mi)^a)
	w[mi < 0] = 0
	metagene = w %*% data
	if(rankBased) metagene = rank(metagene)
	c = 0
	while(c < maxIter){
		mi = getAllMIWz(dataIn, metagene, bin=bin, so=so, negateMI = negateMI)
		delta = sum((mi - premi)^2)
		cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
		print(mi[order(mi, decreasing=T)[1:20]])
		flush.console()
		if(delta < epsilon){
			break
		}
		premi = mi
		mi[mi < 0] = 0
		w = abs(mi)^a / sum(abs(mi)^a)
		w[mi < 0] = 0
		metagene = w %*% data
		if(rankBased) metagene = rank(metagene)
		c = c + 1
	}
	if(c >= maxIter) return (NA)
	return (sort(mi, decreasing=T))
}

findCorrAttractor = function(data, vec, a=10,rankBased = FALSE, maxIter=100, epsilon=1E-14){
	dataIn = data
	if(rankBased){
		vec = rank(vec)
		for(i in 1:nrow(data)){
			dataIn[i,] = rank(data[i,])
		}
	}
	rs = getAllCorWz(dataIn, vec, rankBased=rankBased)
	prers = rs
	w = abs(rs)^a / sum(abs(rs)^a)
	w[rs < 0] = -w[rs < 0]
	metagene = w %*% data
	if(rankBased) metagene = rank(metagene)
	c = 0
	while(c < maxIter){
		rs = getAllCorWz(dataIn, metagene)
		delta = sum((rs - prers)^2)
		cat("Iteration ", (c+1), "\tDelta = ", delta, "\n", sep="")
		print(rs[order(rs, decreasing=T)[1:20]])
		flush.console()
		if(delta < epsilon){
			break
		}
		prers = rs
		w = abs(rs)^a / sum(abs(rs)^a)
		w[rs < 0] = -w[rs < 0]
		metagene = w %*% data
		if(rankBased) metagene = rank(metagene)
		c = c + 1
	}
	if(c >= maxIter) return (NA)
	return (sort(rs, decreasing=T))
}

