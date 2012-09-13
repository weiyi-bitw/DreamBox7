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

CreateMetageneSpace <- function(ge, attractome, map){
  nMeta = length(attractome)
  metaSpace = matrix(NA, nrow=nMeta, ncol=ncol(ge))
  dimnames(metaSpace) = list(names(attractome), colnames(ge))
  pbs = list()
  mappedGenes = map[rownames(ge),1]
  for (i in 1:nMeta){
    #cat(i, "\n")
    #flush.console()
    a = attractome[[i]]
    if(nrow(a) > 10){
      genes = a[1:10,1]
    }else{
      genes = a[, 1]
    }
    il = lapply(genes, function(g){which(mappedGenes == g)})
    ill = sapply(il, length)
    goodIdx = unlist(sapply(il, function(i){ if(length(i) == 1) i}))
    numGood = sum(ill == 1)
    goodMat = NULL
    if(numGood > 0){
      goodMat = ge[goodIdx,]
    }
    badIdx = il[sapply(il, function(i){length(i) > 1 })]

    numBad = length(badIdx)
    
    badMat = NULL
    chosenIdx = NULL
    if(numBad > 0) {
     # geneSum = apply(ge[unlist(il), ],2,sum)
     # chosenIdx = unlist (sapply(badIdx, function(idcs){
	#mis = sapply(idcs, function(i){getMI(geneSum, ge[i,])} )
	#idcs[which.max(mis)]
     # }) )
	#badMat = ge[chosenIdx,]
	badMat = t(sapply(badIdx, function(idcs){
		apply(ge[idcs,], 2, function(x){mean(x, na.rm=TRUE)})
	}) )
    }    
    #pbs[[i]] = c(goodIdx, chosenIdx)
    metaSpace[i,] = (apply(rbind(goodMat, badMat), 2, mean))

    #o = sapply(genes, 
    #          function(g, ge, mappedGenes){
    #           idx = which(mappedGenes == g)
    #          if (length(idx)==1) return (ge[idx,])
    #             if (length(idx)==0) return (rep(NA, ncol(ge)))
    #             return (apply(ge[idx,], 2, function(x) mean(x, na.rm=T)))
    #           },
    #           ge = ge,
    #           mappedGenes = map[rownames(ge), "Gene.Symbol"]
    #           )
    #if(length(genes)==1){metaSpace[i,] = o}
    #else {metaSpace[i,] = apply(o, 1, function(x) mean(x, na.rm=T))}
  }
  #names(pbs) = names(attractome)
  #o = list(metaSpace = metaSpace, pbs = pbs)
  return (metaSpace)
}
