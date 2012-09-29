lazyImputationClnc <- function(c){ # so i'm lazy, bite me!
  for(i in 1:ncol(c)){
    idx = which(is.na(c[,i]))
    if(length(idx) > 0){
      if(class(c[,i])=="numeric") c[idx,i] = mean(c[,i], na.rm=TRUE)
      if(class(c[,i])=="factor") c[idx,i] = 0
    }
  }
  return (c)
}

lazyImputeDFClnc = function(c){
  idx = grep("NOT_IN_OSLOVAL",colnames(c))
  if(length(idx) > 0) c = c[,-idx]
  idx = grep(".Expr", colnames(c))
  if(length(idx) > 0) c = c[,-idx]
  idx = grep("HER2_IHC", colnames(c))
  if(length(idx) > 0) c = c[,-idx]
  for(i in 1:ncol(c)){
    idx = which(is.na(c[,i]))
    if(length(idx) > 0){
      if(class(c[,i]) == "numeric") c[idx, i] = mean(c[-idx,i])
      if(class(c[,i]) == "factor" ){ 
        cc = as.vector(c[,i])
	cc[idx] = "NA"
	levels(c[,i]) <- c(levels(c[,i]), "NA")
	}
    }
  }
  return (c)
}


lazyImputeDFClncFull = function(c){
  c$NOT_IN_OSLOVAL_stage = factor(c$NOT_IN_OSLOVAL_stage)
  idx = grep(".Expr", colnames(c))
  if(length(idx) > 0) c = c[,-idx]
  killIdx = NULL
  for(i in 1:ncol(c)){
    idx = which(is.na(c[,i]))
    if(length(idx) > 0.3 * nrow(c)) killIdx = c(killIdx, i)
  }
  c = c[, -killIdx]
  for(i in 1:ncol(c)){
    idx = which(is.na(c[,i]))
    if(length(idx) > 0){
      if(class(c[,i]) == "numeric") c[idx, i] = mean(c[-idx,i])
      if(class(c[,i]) == "factor" ){
        cc = as.vector(c[,i])
        cc[idx] = "NA"
        c[, i] = factor(cc)
        }
    }
  }
  return (c)

}

CreateMetageneSpace <- function(ge, attractome, map, chosenProbes = NULL){
  if(is.null(chosenProbes)) {
  nMeta = length(attractome)
  metaSpace = matrix(NA, nrow=nMeta, ncol=ncol(ge))
  dimnames(metaSpace) = list(names(attractome), colnames(ge))
  pbs = list()
  mappedGenes = rep(NA, nrow(ge))
  names(mappedGenes) = rownames(ge)
  mappedGenes[rownames(map)] = map[,1]
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
    goodIdx = sapply(il, function(i){ if(length(i) == 1) i})
    goodIdx = goodIdx[sapply(goodIdx, function(x){!is.null(x)})]
    numGood = sum(ill == 1)
    goodMat = NULL
    if(numGood > 0){
      goodMat = ge[unlist(goodIdx),]
    }
    badIdx = il[sapply(il, function(i){length(i) > 1 })]

    numBad = length(badIdx)
    
    badMat = NULL
    chosenIdx = NULL
    if(numBad > 0) {
      geneSum = apply(ge[unlist(il), ],2,sum)
      chosenIdx = lapply(badIdx, function(idcs){
	mis = sapply(idcs, function(i){getMI(geneSum, ge[i,])} )
	idcs[which(mis > 0.5)]
      })
      chosenIdx = chosenIdx[sapply(chosenIdx, function(x){length(x)>0})]
	#badMat = ge[chosenIdx,]
	badMat = t(sapply(chosenIdx, function(idcs){
		if(length(idcs) > 1){
			apply(ge[idcs,], 2, function(x){mean(x, na.rm=TRUE)})
		}else if(length(idcs) == 1){
			ge[idcs,]
		}else{
			rep(NA, ncol(ge))
		}
	}) )
	if(length(chosenIdx) == 0) {chosenIdx = NULL; badMat = NULL}
    }    
    pbs[[i]] = c(goodIdx, chosenIdx)
    metaSpace[i,] = (apply(rbind(goodMat, badMat), 2, function(x){mean(x, na.rm=TRUE)}))

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
  names(pbs) = names(attractome)
  o = list(metaSpace = metaSpace, pbs = pbs)
  return (o)
  }else{

	metaSpace = t(sapply(chosenProbes, function(pb){
		gmat = sapply(pb, function(p, ge){
			if(length(p) > 1){
				apply(ge[p,], 2, mean)
			}else{
				ge[p,]
			}
		}, ge = ge)
		if(!is.null(dim(gmat))) {apply(gmat, 1, mean)}
		else{gmat}
	}) )
	return(metaSpace)

  }
}
CreateGeneSpace <- function(ge, oncogenes, map){
  ng = length(oncogenes)
  gSpace = matrix(NA, nrow=ng, ncol=ncol(ge))
  dimnames(gSpace) = list(oncogenes, colnames(ge))
  mappedGenes = map[rownames(ge), "Gene.Symbol"]
  for(i in 1:ng){
    g = oncogenes[i]
    idx = which(mappedGenes == g)
    if(length(idx)==1){gSpace[i,] = ge[idx,]}
    else{gSpace[i,] = apply(ge[idx,], 2, mean)}
  }
  return (gSpace)
}

expandClnc = function(c){
  h.IDC =as.numeric(c$histological_type=="IDC")
  h.ILC =as.numeric(c$histological_type=="ILC")
  h.IDCpILC =as.numeric(c$histological_type=="IDC+ILC")
  h.IDCnMED =as.numeric(c$histological_type=="IDC-MED")
  h.IDCnMUC =as.numeric(c$histological_type=="IDC-MUC")
  h.IDCnTUB =as.numeric(c$histological_type=="IDC-TUB")
  #h.MIXED =as.numeric(c$histological_type=="MIXED NST AND A SPECIAL TYPE")
  #h.OTHER =as.numeric(c$histological_type=="OTHER")
  #h.OTHERINV =as.numeric(c$histological_type=="OTHER INVASIVE")
  #h.INVTUMOR =as.numeric(c$histological_type=="INVASIVE TUMOR")
  h.other = as.numeric(c$histological_type=="MIXED NST AND A SPECIAL TYPE" | c$histological_type=="OTHER" | c$histological_type=="OTHER INVASIVE" | c$histological_type=="INVASIVE TUMOR" | c$histological_type=="PHYL")

  er.P=as.numeric(c$ER_IHC_status=="pos")
  er.N=as.numeric(c$ER_IHC_status=="neg")
  
  tr.CT = as.numeric((c$Treatment == "CT") | (c$Treatment == "CT/HT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT"))
  tr.HT = as.numeric((c$Treatment == "HT") | (c$Treatment == "CT/HT") | (c$Treatment == "HT/RT") | (c$Treatment =="CT/HT/RT"))
  tr.RT = as.numeric((c$Treatment == "RT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT") | (c$Treatment == "HT/RT"))
  trMat = cbind(tr.CT, tr.HT, tr.RT)

  gd.1 = as.numeric(c$grade==1)
  gd.2 = as.numeric(c$grade==2)
  gd.3 = as.numeric(c$grade==3)

  #her2.1 = as.numeric(c$HER2_IHC_status==1)
  #her2.2 = as.numeric(c$HER2_IHC_status==2)
  #her2.3 = as.numeric(c$HER2_IHC_status==3)
  
  her2.snp6.gain = as.numeric(c$HER2_SNP6_state=="GAIN")
  her2.snp6.loss = as.numeric(c$HER2_SNP6_state=="LOSS")
  her2.snp6.neut = as.numeric(c$HER2_SNP6_state=="NEUT")
  
  cmat<-data.frame(c[, c(1:3)], gd.1, gd.2, gd.3, h.IDC, h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, h.other,# h.MIXED, h.OTHER, h.OTHERINV, h.INVTUMOR,
                   er.N, er.P, tr.CT, tr.HT, tr.RT, her2.snp6.gain, her2.snp6.loss, her2.snp6.neut)#, p.LumA, p.LumB, p.Her2, p.Normal, p.Basal, her2.1, her2.2, her2.3)
  for(i in 4:ncol(cmat)){
    cmat[,i] = factor(cmat[,i])
  }
  return(cmat)
}


expandClncFull = function(c){
  h.IDC =as.numeric(c$histological_type=="IDC")
  h.ILC =as.numeric(c$histological_type=="ILC")
  h.IDCpILC =as.numeric(c$histological_type=="IDC+ILC")
  h.IDCnMED =as.numeric(c$histological_type=="IDC-MED")
  h.IDCnMUC =as.numeric(c$histological_type=="IDC-MUC")
  h.IDCnTUB =as.numeric(c$histological_type=="IDC-TUB")
  #h.MIXED =as.numeric(c$histological_type=="MIXED NST AND A SPECIAL TYPE")
  #h.OTHER =as.numeric(c$histological_type=="OTHER")
  #h.OTHERINV =as.numeric(c$histological_type=="OTHER INVASIVE")
  #h.INVTUMOR =as.numeric(c$histological_type=="INVASIVE TUMOR")
  h.other = as.numeric(c$histological_type=="MIXED NST AND A SPECIAL TYPE" | c$histological_type=="OTHER" | c$histological_type=="OTHER INVASIVE" | c$histological_type=="INVASIVE TUMOR" | c$histological_type=="PHYL")

  er.P=as.numeric(c$ER_IHC_status=="pos")
  er.N=as.numeric(c$ER_IHC_status=="neg")
  
  tr.CT = as.numeric((c$Treatment == "CT") | (c$Treatment == "CT/HT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT"))
  tr.HT = as.numeric((c$Treatment == "HT") | (c$Treatment == "CT/HT") | (c$Treatment == "HT/RT") | (c$Treatment =="CT/HT/RT"))
  tr.RT = as.numeric((c$Treatment == "RT") | (c$Treatment == "CT/HT/RT") | (c$Treatment == "CT/RT") | (c$Treatment == "HT/RT"))
  trMat = cbind(tr.CT, tr.HT, tr.RT)

  gd.1 = as.numeric(c$grade==1)
  gd.2 = as.numeric(c$grade==2)
  gd.3 = as.numeric(c$grade==3)

  grp.1 = as.numeric(c$NOT_IN_OSLOVAL_group==1)
  grp.2 = as.numeric(c$NOT_IN_OSLOVAL_group==2)
  grp.3 = as.numeric(c$NOT_IN_OSLOVAL_group==3)
  grp.4 = as.numeric(c$NOT_IN_OSLOVAL_group==4)

  stg.0 = as.numeric(c$NOT_IN_OSLOVAL_stage==0)
  stg.1 = as.numeric(c$NOT_IN_OSLOVAL_stage==1)
  stg.2 = as.numeric(c$NOT_IN_OSLOVAL_stage==2)
  stg.3 = as.numeric(c$NOT_IN_OSLOVAL_stage==3)
  stg.4 = as.numeric(c$NOT_IN_OSLOVAL_stage==4)

  #her2.1 = as.numeric(c$HER2_IHC_status==1)
  #her2.2 = as.numeric(c$HER2_IHC_status==2)
  #her2.3 = as.numeric(c$HER2_IHC_status==3)
  
  her2.snp6.gain = as.numeric(c$HER2_SNP6_state=="GAIN")
  her2.snp6.loss = as.numeric(c$HER2_SNP6_state=="LOSS")
  her2.snp6.neut = as.numeric(c$HER2_SNP6_state=="NEUT")
  
  st.1 = as.numeric(c$NOT_IN_OSLOVAL_Site==1)
  st.2 = as.numeric(c$NOT_IN_OSLOVAL_Site==2)
  st.3 = as.numeric(c$NOT_IN_OSLOVAL_Site==3)
  st.4 = as.numeric(c$NOT_IN_OSLOVAL_Site==4)
  st.5 = as.numeric(c$NOT_IN_OSLOVAL_Site==5)


  cmat<-data.frame(c[, c(1:3)], gd.1, gd.2, gd.3, h.IDC, h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, h.other,# h.MIXED, h.OTHER, h.OTHERINV, h.INVTUMOR,
                   er.N, er.P, tr.CT, tr.HT, tr.RT, her2.snp6.gain, her2.snp6.loss, her2.snp6.neut, grp.1, grp.2, grp.3, grp.4, stg.0, stg.1, stg.2, stg.3, stg.4, st.1, st.2, st.3, st.4, st.5)#, p.LumA, p.LumB, p.Her2, p.Normal, p.Basal, her2.1, her2.2, her2.3)
  for(i in 4:ncol(cmat)){
    cmat[,i] = factor(cmat[,i])
  }
  return(cmat)
}
