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
	c[, i] = factor(cc)
	}
    }
  }
  return (c)
}
lazyImputeDFClncFULL = function(c){
  c$NOT_IN_OSLOVAL_stage <- factor(c$NOT_IN_OSLOVAL_stage)
  class(c$HER2_IHC_status) = "factor"
  killIdx = NULL
  for(i in 1:ncol(c)){
    mis = sum(is.na(c[,i]))
    if(mis > (nrow(c) * 0.3) ) killIdx = c(killIdx , i)
  }
  c = c[,-killIdx]
  idx = grep(".Expr", colnames(c))
  if(length(idx) > 0) c = c[,-idx]

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
  idx = intersect(rownames(map) , rownames(ge))
  mappedGenes[idx] = map[idx,1]
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
    pbs[[i]] = sapply(c(goodIdx, chosenIdx), names)
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
		pb = sapply(pb, function(p){intersect(p, rownames(ge))})
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
  mappedGenes = rep(NA, nrow(ge))
  names(mappedGenes) = rownames(ge)
  idx = intersect(rownames(ge) , rownames(map))
  mappedGenes[idx] = map[idx,1]
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

expandClncFULL = function(c){
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
  stg.3 = as.numeric(c$NOT_IN_OSLOVAL_stage==3 | c$NOT_IN_OSLOVAL_stage==4)

#  her2.1 = as.numeric(c$HER2_IHC_status==1)
#  her2.2 = as.numeric(c$HER2_IHC_status==2)
#  her2.3 = as.numeric(c$HER2_IHC_status==3)
  
  her2.snp6.gain = as.numeric(c$HER2_SNP6_state=="GAIN")
  her2.snp6.loss = as.numeric(c$HER2_SNP6_state=="LOSS")
  her2.snp6.neut = as.numeric(c$HER2_SNP6_state=="NEUT")
  
  st.1 = as.numeric(c$NOT_IN_OSLOVAL_Site==1)
  st.2 = as.numeric(c$NOT_IN_OSLOVAL_Site==2)
  st.3 = as.numeric(c$NOT_IN_OSLOVAL_Site==3)
  st.4 = as.numeric(c$NOT_IN_OSLOVAL_Site==4)
  st.5 = as.numeric(c$NOT_IN_OSLOVAL_Site==5)

  mnpsl.pre = as.numeric(c$NOT_IN_OSLOVAL_menopausal_status_inferred=="pre" )
  mnpsl.post = as.numeric(c$NOT_IN_OSLOVAL_menopausal_status_inferred=="post" )

#  cell.hi = as.numeric(c$NOT_IN_OSLOVAL_cellularity == "high")
#  cell.mod = as.numeric(c$NOT_IN_OSLOVAL_cellularity == "moderate")
#  cell.lo = as.numeric(c$NOT_IN_OSLOVAL_cellularity == "low")

  p.LumA = as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype == "LumA")
  p.LumB = as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype == "LumB")
  p.Basal = as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype == "Basal")
  p.Her2 = as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype == "Her2")
  p.Normal = as.numeric(c$NOT_IN_OSLOVAL_Pam50Subtype == "Normal")

#  p53.mut = as.numeric(c$NOT_IN_OSLOVAL_P53_mutation_status == "MUT")
#  p53.wt = as.numeric(c$NOT_IN_OSLOVAL_P53_mutation_status == "WT")

#  ic.1 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==1)
#  ic.2 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==2)
#  ic.3 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==3)
#  ic.4 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==4)
#  ic.5 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==5)
#  ic.6 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==6)
#  ic.7 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==7)
#  ic.8 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==8)
#  ic.9 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==9)
#  ic.10 = as.numeric(c$NOT_IN_OSLOVAL_IntClustMemb==10)

#  gf.nn = as.numeric(c$NOT_IN_OSLOVAL_Genefu == "ER-/HER2-")
#  gf.pnhi = as.numeric(c$NOT_IN_OSLOVAL_Genefu == "ER+/HER2- High Prolif")
#  gf.pnlo = as.numeric(c$NOT_IN_OSLOVAL_Genefu == "ER+/HER2- Low Prolif")
#  gf.p = as.numeric(c$NOT_IN_OSLOVAL_Genefu == "HER2+")


  cmat<-data.frame(c[, c(1:3, 12:13)], gd.1, gd.2, gd.3, h.IDC, h.ILC, h.IDCpILC, h.IDCnMED, h.IDCnMUC, h.IDCnTUB, h.other,# h.MIXED, h.OTHER, h.OTHERINV, h.INVTUMOR,
                   er.N, er.P, tr.CT, tr.HT, tr.RT, her2.snp6.gain, her2.snp6.loss, her2.snp6.neut, 
		grp.1, grp.2, grp.3, grp.4, stg.0, stg.1, stg.2, stg.3, st.1, st.2, st.3, st.4, st.5, p.LumA, p.Basal, p.LumB, p.Normal, p.Her2,
		mnpsl.pre, mnpsl.post)
		#, cell.hi, cell.mod, cell.lo)
		#, p53.mut, p53.wt, her2.1, her2.2, her2.3)
		#, gf.nn, gf.pnhi, gf.pnlo, gf.p)
		#, ic.1, ic.2, ic.3, ic.4, ic.5, ic.6, ic.7, ic.8, ic.9, ic.10)
  for(i in 6:ncol(cmat)){
    cmat[,i] = factor(cmat[,i])
  }
  return(cmat)
}

getTags = function(ft){
	taglist = strsplit(ft, "\\.")
	tags = sapply(taglist, function(x){x[1]})
	return (tags)
}

removeTaggedFeatures = function(colName, ft){
	t = getTags(colName)
	tags = getTags(ft)
	killIdx = which(t %in% tags)
	out = colName[-killIdx]
	return (out)
}

