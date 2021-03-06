# File name: call_c.r

# dyn.load("call_c.so")
# For Windows will like this
# dyn.load("C:/Windows/Desktop/call_c.dll")

# symbol.C("hello")
# is.loaded(symbol.C("hello"))

# a <- 1 : 9
getMI = function(x, y, bin=6, so=3, rankBased=FALSE, normalize=TRUE, negateMI = FALSE){
  n <- length(x)
  if(length(y) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  if(is.factor(y)){
    yph = y
    idx = which(!is.na(yph))
    x = x[idx]
    yph = yph[idx]

    y = as.numeric(yph)
    l = levels(yph)
    nl = length(l)
    if(is.ordered(yph)){soph = 3}
    else {soph = 1}
    out = .C("mi2DiffBins", x = as.double(x), y = as.double(y), n = as.integer(length(yph)), binx = as.integer(bin), biny = as.integer(nl), sox = as.integer(so), soy = as.integer(soph), miOut = 0, norm = as.integer(normalize), negateMI = as.integer(negateMI))
  }else{

    if(rankBased){
      x = rank(x)
      y = rank(y)
    }
      out = .C("mi2", x = as.double(x), y=as.double(y), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), miOut = 0, norm = as.integer(normalize), negateMI = as.integer(negateMI))

  }

  mi = out$miOut
  return (mi)

}
getMI3 = function(x, y, z, bin=6, so=3, rankBased=FALSE){
  n <- length(x)
  if(length(y) != n | length(z) != n){stop("legnth of three vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
    if(rankBased){
      x = rank(x)
      y = rank(y)
      z = rank(z)
    }
      out = .C("mi3", x = as.double(x), y=as.double(y), z=as.double(z), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), miOut = 0)

  mi = out$miOut
  return (mi)

}



getMI2vs1 = function(x, y, z, bin=6, so=3, normalize=FALSE){
  n <- length(x)
  if(length(y) != n | length(z) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}
  
  out = .C("mi2vs1", x = as.double(x), y=as.double(y), z=as.double(z), n=as.integer(n), bin=as.integer(bin), so = as.integer(so), mi = 0, norm=as.integer(normalize))

  mi = out$mi
  return (mi)
  
}

getCCDIdx = function(predictions, observations){
	out = .Call("ccdiR2C", predictions, observations)
	return (out); 
}

getWCCDIdx = function(predictions, w, observations){
	out = .Call("ccdiwR2C", predictions, w, observations)
	return (out); 
}

getAllCCDIWz = function(x, observations, sorted=FALSE){
        n = ncol(x)
        m = nrow(x)
        out = .C("concordance_index_all", predictions=as.double(x), observations=as.double(observations), mIn = as.integer(m), nIn = as.integer(n), score = rep(0, m))
        score = out$score
        names(score) = rownames(x)
        if(sorted) score = sort(score, decreasing=T)
        return (score)
}


getCorr = function(x, y, rankBased = FALSE){
	n = length(x)
	if(length(y) != n){stop("length of two vectors are different!")}
	if(rankBased){
		x = rank(x)
		y = rank(y)
	}

	out = .C("cor", x=as.double(x), y=as.double(y), n=as.integer(n), cOut=as.double(0))

	r = out$cOut
	return(r)
}

getAllMIWz = function(x, vec, bin=6, so=3, rankBased=FALSE, normalize = TRUE, sorting = FALSE, negateMI = FALSE){
  m = nrow(x)
  n = ncol(x)
  if(length(vec) != n){stop("legnth of two vectors are different!")}
  if(so >= bin){stop("spline order must be less than bin")}

  if(rankBased){
    for(i in 1:m){
      x[i,] = rank(x[i,])
    } 
    vec = rank(vec)
  }
  garbage = rep(-1, m)
  out = .C("getAllMIWz", data = as.double(x), vec = as.double(vec), mi = as.double(garbage), m = as.integer(m), n = as.integer(n), bin = as.integer(bin),so = as.integer(so), norm=as.integer(normalize), negateMI = as.integer(negateMI))

  mi = out$mi
  names(mi) = rownames(x)
  if(sorting) mi = sort(mi, decreasing=TRUE)
  return (mi)
}

getAllPWMIWz = function(x, bin=6, so=3, rankBased=FALSE, normalize = TRUE, negateMI = FALSE){
  m = nrow(x)
  n = ncol(x)
  if(so >= bin){stop("spline order must be less than bin")}

  if(rankBased){
    for(i in 1:m){
      x[i,] = rank(x[i,])
    } 
    vec = rank(vec)
  }
  garbage = rep(-1, m*(m-1)/2)
  out = .C("getAllPWMIWz", data = as.double(x), mi = as.double(garbage), m = as.integer(m), n = as.integer(n), bin = as.integer(bin),so = as.integer(so), norm=as.integer(normalize), negateMI = as.integer(negateMI))

  mi = out$mi
  return (mi)
}

getAllCorWz = function(x, vec, rankBased=FALSE, sorting=FALSE){
  m = nrow(x)
  n = ncol(x)
  if(length(vec) != n){stop("legnth of two vectors are different!")}

  if(rankBased){
    for(i in 1:m){
      x[i,] = rank(x[i,])
    } 
    vec = rank(vec)
  }
  garbage = rep(-999, m)
  out = .C("getAllCorWz", data = as.double(x), vec = as.double(vec), m = as.integer(m), n = as.integer(n), rs=as.double(garbage))

  rs = out$rs
  names(rs) = rownames(x)
  if(sorting) rs = sort(rs, decreasing=TRUE)
  return (rs)
}

logRankScore = function(survobj){
	if(class(survobj) != "Surv"){stop("Input survobj must be a Surv object!!")}
	t = survobj[,1]
	e = survobj[,2]
	n = length(t)

	o = rep(0, n)
	lrScore = rep(0, n)

	out = .C("logRankScore", x=as.double(t), evt=as.integer(e), nIn=as.integer(n), o=as.integer(o), lrScore=as.double(lrScore))
	return (out$lrScore)
}

testNA = function(i){
	out = .Call("testNA", i)
	return (out)
}

lowerTriIndex = function(m, diag=FALSE){
	if(diag){
		n = m * (m+1) / 2
	}else{
		n = m * (m-1) / 2
	}
	x = rep(-1, n)
	y = rep(-1, n)
	out = .C("lowerTriIndex", x=as.integer(x), y=as.integer(y), m = as.integer(m), diag=as.integer(diag))

	return(list(x=out$x+1, y=out$y+1))	
}

weihouette = function(mis, cluster.member){
	k = length(table(cluster.member))
	m = length(cluster.member)
	lenpw = length(mis)
	w = rep(0, k)
	a = rep(0, k)
	b = rep(0, k)

	out = .C("weihouette", mis = as.double(mis), clust_member = as.integer(cluster.member-1), k = as.integer(k), m = as.integer(m), lenpw = as.integer(lenpw), w = as.double(w), a = as.double(a), b = as.double(b) )

	return (list(w = out$w, a = out$a, b = out$b))

}

# test.c(a)

# dyn.unload("call_c.so")
# For Windows will like this
# dyn.unload("C:/Windows/Desktop/call_c.dll")

