getTumorCounts <- function(vcf){
	sapply(grep("(F|R).Z", names(geno(vcf)), value=TRUE), function(i) geno(vcf)[[i]][,"TUMOUR"])
}

#' Extract tumour depth
#' @param vcf 
#' @return numeric()
#' 
#' @author mg14
#' @export
getTumorDepth <- function(vcf){
	if("t_alt_count" %in% colnames(info(vcf))){ ## consensus data, snv and indel
		info(vcf)$t_alt_count + info(vcf)$t_ref_count
	}else{ ## older data
		if("FAZ" %in% rownames(geno(header(vcf)))){
			rowSums(getTumorCounts(vcf))
		}else{
			geno(vcf)$DEP[,2]
		}
	}
}


#' Get alt alleles counts 
#' @param vcf 
#' @return numeric()
#' 
#' @author mg14
#' @export
getAltCount <- function(vcf){
	if("t_alt_count" %in% colnames(info(vcf))){ ## consensus data, snv and indel
		info(vcf)$t_alt_count
	}else{ ## older formats
		if("FAZ" %in% rownames(geno(header(vcf)))){ ## ie subs
			c <- getTumorCounts(vcf)
			t <- c[,1:4] + c[,5:8]
			colnames(t) <- substring(colnames(t),2,2)
			a <- as.character(unlist(alt(vcf)))
			a[!a%in%c('A','T','C','G')] <- NA
			sapply(seq_along(a), function(i) if(is.na(a[i])) NA else t[i, a[i]])
		}
		else{ ## ie indel
			#(geno(vcf)$PP + geno(vcf)$NP + geno(vcf)$PB + geno(vcf)$NB)[,"TUMOUR"]
			geno(vcf)$MTR[,2]
		}
	}
}

#' @importFrom VGAM pbetabinom
dtrbinom <- function(x, size, prob, xmin=0) dbinom(x,size,prob) / pbinom(xmin-1, size, prob, lower.tail=FALSE)
pbetabinom <- function(x, size, prob, rho){
	if(rho!=0)
		VGAM::pbetabinom(x, size, prob, rho)
	else
		pbinom(x, size, prob)
}
#' @importFrom VGAM dbetabinom
dbetabinom <- function(x, size, prob, rho){
	if(rho!=0)
		VGAM::dbetabinom(x, size, prob, rho)
	else
		dbinom(x, size, prob)
}
dtrbetabinom <- function(x, size, prob, rho, xmin=0) {y <- dbetabinom(x,size,prob,rho) / (1-pbetabinom(xmin-1, size, prob, rho))
	y[x<xmin] <- 0
	return(y)}
ptrbetabinom <- function(x, size, prob, rho, xmin=0) {
	pmin <- pbetabinom(xmin-1, size, prob, rho)
	pmax(0,(pbetabinom(x,size,prob,rho) - pmin) / (1-pmin))}


mergeClusters <- function(clusters, deltaFreq=0.05){
	if(nrow(clusters) <= 1) return(clusters)
	h <- hclust(dist(clusters$proportion))
	ct <- cutree(h, h=deltaFreq)
	cp <- as.matrix(cophenetic(h))
	Reduce("rbind",lapply(unique(ct), function(i) {
						n_ssms <- sum(clusters$n_ssms[ct==i])
						w <- max(cp[ct==i,ct==i])
						data.frame(new.cluster=i, n_ssms=n_ssms, proportion=sum(clusters$proportion[ct==i]*clusters$n_ssms[ct==i])/n_ssms, width=w)
					}
			))
}

#' Compute possible Mutation copy number states
#' @param cn Copy number with consensus annotation meta data
#' @param clusters 
#' @param purity 
#' @param gender 
#' @param isWgd 
#' @param deltaFreq
#' @return list of length nrow(bb), can be added to mcols(bb)
#' 
#' @author mg14
#' @export
defineMcnStates <- function(cn, clusters, purity, gender='female', isWgd= FALSE, deltaFreq=0.05){
	P <- vector(mode='list', length(cn))
	uniqueBB <- unique(cn)
	overlaps <- findOverlaps(uniqueBB, cn)
	
	majorCN <- split(cn$major_cn[subjectHits(overlaps)], queryHits(overlaps))
	m <- cn$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
	if(gender=='male')
		m[as.character(seqnames(cn)) %in% c('X','Y')] <- 0
	minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps))	
	h <- selectHits(overlaps, "first")
	H <- selectHits(overlaps, "last")
	
	cnNormal <- 2 - (gender=='male' & seqnames(uniqueBB)=="X" | seqnames(uniqueBB)=="Y")
	
	# Fix cluster and purity discrepancies
	clusters$proportion[which.max(clusters$proportion)] <- purity
	
	cloneFreq <- split(cn$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
	cnStates <- matrix(0, nrow=10000, ncol=6)
	colnames(cnStates) <- c("state","m","f","n.m.s","pi.m.s","s")
	
	power.c <- rep(0, nrow(clusters))
	
	
	for( i in seq_along(uniqueBB)){
		if(!i %in% names(majorCN)) next
		majcni <- majorCN[[as.character(i)]]
		mincni <- minorCN[[as.character(i)]]
		cfi <- cloneFreq[[as.character(i)]]
		effCnTumor <- sum((majcni + mincni)*cfi)
		effCnNormal <- as.numeric(cnNormal[i]) * (1-purity)
		
		if(any(is.na(majcni))) next
		
		mixFlag <- FALSE
		subclonalGainFlag <- FALSE
		clonalFlag <- TRUE
		majdelta <- 0
		mindelta <- 0
		
		majanc <- majder <- majcni
		minanc <- minder <- mincni
		
		if(length(cfi)>1){ # multiple (subclonal) CN states, if so add clonal option (ie. mixture of both states), subclonal states only change by 1..delta(CN)
			d <- colSums(abs(rbind(majcni, mincni) - c(1,1) * (1+ isWgd)))
			derived <- d == max(d) ## derived state further from diploid/tetraploid
			if(all(derived)) next 
			majanc <- majcni[!derived]
			minanc <- mincni[!derived]
			majder <- majcni[derived]
			minder <- mincni[derived]
			majdelta <- majcni[derived] - majcni[!derived]
			mindelta <- mincni[derived] - mincni[!derived]
			majcni <- c(min(majcni), # clonal, sub on allele that doesn't change
					(majcni[!derived] + cfi[derived]/purity*majdelta), # clonal, sub on allele that does change
					(majcni[derived] >0) + (majdelta > 0)) # subclonal, subs after or before CNA, m=1,delta
			mincni <- c(min(mincni), (mincni[!derived] + cfi[derived]/purity*mindelta), (mincni[derived] >0) + (mindelta > 0))
			majdelta <- c(0, cfi[derived]/purity*majdelta, majdelta)
			mindelta <- c(0, cfi[derived]/purity*mindelta, mindelta)
			cfi <- c(purity, purity,  cfi[derived])
			mixFlag <- c(FALSE, TRUE, FALSE)
			clonalFlag <- c(TRUE,TRUE,FALSE)
			subclonalGainFlag <- c(FALSE, FALSE, TRUE)
		}
		
		
		a <- sapply(clusters$proportion, function(p) all(abs(p-cfi) > deltaFreq)) # subclone(s) not coinciding with CN change
		if(any(a)){ # assume subclones have derived from most abundant CN state
			majcni <- c(majcni, rep(majcni[which.max(cfi)]>0, sum(a))+0)
			mincni <- c(mincni, rep(mincni[which.max(cfi)]>0, sum(a))+0)
			cfi <- c(cfi, clusters$proportion[a])
			mixFlag <- c(mixFlag, rep(FALSE, sum(a)))
			clonalFlag <- c(clonalFlag, rep(FALSE, sum(a)))
			subclonalGainFlag <- c(subclonalGainFlag, rep(FALSE, sum(a)))
			majdelta <- c(majdelta, rep(0,sum(a)))
			mindelta <- c(mindelta, rep(0,sum(a)))
		}
		totcni <- majcni+mincni
		if(all(totcni==0)) next
		
		k <- 0
		for( j in seq_along(majcni)){
			if(majcni[j]==0 & mincni[j]==0) {
				f <- m <- 0 # allele frequency
				pi.m.s <- n.m.s <- 1
				l <- 1
			}else{
				l <- 1:max(majcni[j], mincni[j]) # mincni>majcni can occur if minor allele changes subclonally
				m <- l
				n.m.s <- rep(1, length(l)) #multiplicity of cn states
				if(!mixFlag[j]){ # single subclone, or no subclonal cn
					f <- m * cfi[j] / (effCnTumor + effCnNormal)
					if(mincni[j] > 0)
						n.m.s[1:min(majcni[j], mincni[j])] <- 2
					pi.m.s <- n.m.s/sum(n.m.s)
				}else{ # coexisting cn subclones, use mixture
					d <- if(majdelta[j] != 0) majdelta[j] else mindelta[j] 
					M <- max(mincni[j]*(mindelta[j]!=0), majcni[j]*(majdelta[j]!=0))
					m <- (d > 0):M + (M-floor(M))
					l <- seq_along(m)
					f <- m *cfi[j] / (effCnTumor + effCnNormal) # Variant on major allele
					n.m.s <- rep(1, length(l))
					pi.m.s <- n.m.s/sum(n.m.s)
				}
			}
			cnStates[k + l,"state"]=j
			cnStates[k + l,"m"]=m
			cnStates[k + l,"f"]=f
			cnStates[k + l,"pi.m.s"]=pi.m.s
			cnStates[k + l,"n.m.s"]=n.m.s
			k <- k + length(l)
		}
		whichStates <- (1:k)[cnStates[1:k,"f"]>0]
		
		# State probabilities - based on cell fractions
		cfi.s <- unique(cfi)
		pi.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
		pi.s <- pi.s/sum(pi.s)
		
		c.to.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, which.min(abs(clusters$proportion - p)), NA)) # map to cluster
		s.to.c <- sapply(clusters$proportion, function(c) which.min(abs(cfi.s - c)))
		
		cnStates[1:k,"s"] = as.numeric(factor(cfi, levels=cfi.s))[cnStates[1:k,"state"]]
		
		timing_param <- cbind(cnStates[whichStates,,drop=FALSE], cfi=cfi[cnStates[whichStates,"state"]], pi.s=pi.s[cnStates[whichStates,"s"]], P.m.sX=NA,P.m.sX.lo=NA, P.m.sX.up=NA, T.m.sX=NA, T.m.sX.lo=NA, T.m.sX.up=NA, power.s=NA, power.m.s = NA,
				majCN=majcni[cnStates[whichStates,"state"]], minCN=mincni[cnStates[whichStates,"state"]], 
				majDelta = majdelta[cnStates[whichStates,"state"]], minDelta = mindelta[cnStates[whichStates,"state"]], 
				clonalFlag=clonalFlag[cnStates[whichStates,"state"]], subclonalGainFlag=subclonalGainFlag[cnStates[whichStates,"state"]], mixFlag=mixFlag[cnStates[whichStates,"state"]], majCNanc=majanc, minCNanc=minanc, majCNder=majder, minCNder=minder)
		P[[h[i]]] <- timing_param
		if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]
		
	}
	return(P)
} 


#' Compute timing parameters
#' @param vcf A vcf object of ssnms. See VariantAnnotation::readVcf()
#' @param bb The copy number as a GRanges() object, meta data in consensus format. See loadBB()
#' @param clusters A data.frame with the cluster proportion and n_ssms
#' @param purity The purity of the samples
#' @param gender 'male' or 'female'
#' @param isWgd TRUE/FALSE 
#' @param xmin min read support. Needed for power calculations
#' @param rho Dispersion parameter
#' @param n.boot Number of bootstrap samples for confidence intervals.
#' @param deltaFreq The difference between subclonal SNV and CN states within which they'd be merged
#' @return A list with elements (D: Annotated Data.Frame, can be added to vcf object; P: Timing parameters to be added to CN Ranges; power.c power of each cluster).
#' @example inst/example/example.R
#' 
#' @author mg14
#' @export
computeMutCn <- function(vcf, bb, clusters=data.frame(cluster=1, proportion=max(bb$clonal_frequency,na.rm=TRUE), n_ssms=100), purity=max(bb$clonal_frequency,na.rm=TRUE), gender='female', isWgd= FALSE, xmin=3, rho=0, n.boot=200, deltaFreq=0.05){
	n <- nrow(vcf)
	D <- DataFrame(MutCN=rep(NA,n), MutDeltaCN=rep(NA,n), MajCN=rep(NA,n), MinCN=rep(NA,n), MajDerCN=rep(NA,n), MinDerCN=rep(NA,n), CNF=rep(NA,n), CNID =as(vector("list", n),"List"), pMutCN=rep(NA,n), pGain=rep(NA,n),pSingle=rep(NA,n),pSub=rep(NA,n), pAllSubclones=as(vector(mode="list",n),"List"), pMutCNTail=rep(NA,n))	
	P <- defineMcnStates(bb,clusters, purity, gender, isWgd, deltaFreq=deltaFreq)
	if(n==0)
		return(list(D=D, P=P))
	
	altCount <- getAltCount(vcf)
	tumDepth <- getTumorDepth(vcf)
	names(altCount) <- names(tumDepth) <- NULL
	
	# Match VCF and CN
	overlaps <- findOverlaps(vcf, bb)
	D[["CNID"]] <- as(overlaps, "List")
	majorCN <- split(bb$major_cn[subjectHits(overlaps)], queryHits(overlaps))
	m <- bb$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
	if(gender=='male')
		m[as.character(seqnames(bb)) %in% c('X','Y')] <- 0
	minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps))	
	h <- selectHits(overlaps, "first")
	H <- selectHits(overlaps, "last")
	
	cnNormal <- 2 - (gender=='male' & seqnames(vcf)=="X" | seqnames(vcf)=="Y")
	
	# Fix cluster and purity discrepancies
	clusters$proportion[which.max(clusters$proportion)] <- purity
	
	cloneFreq <- split(bb$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
	
	power.c <- rep(0, nrow(clusters))
	
	deltaFreq <- 0.05 # merge clusters withing deltaFreq
	
	boundaryHits <- countOverlaps(vcf, unique(bb)) > 1 # indels overlapping with CN boundaries
	
	for(globalIt in 1:2){ # 2 iterations, fist local (ie per segment) fit, then global
		for( i in which( (diff(c(-1, h)) !=0 | is.na(diff(c(-1, h)) !=0) ) & ! boundaryHits )){
			if(!i %in% names(majorCN)) next
			if(is.na(h[i])) next
			if(if(i>1) h[i] != h[i-1] | is.na(h[i-1]) else TRUE){ #ie. new segment
				
				hh <- which(h==h[i] & !is.na(altCount) &! is.na(tumDepth))
				if(length(hh)==0) next
				
				if(is.null(bb$timing_param[[h[i]]])){
					cnStates <- P[[h[i]]]
					if(is.null(cnStates)) next
					whichStates <- 1:nrow(cnStates)
					
					# State probabilities - based on cell fractions
					cfi.s <- unique(cnStates[,"cfi"])
					pi.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
					pi.s <- pi.s/sum(pi.s)
					
					c.to.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, which.min(abs(clusters$proportion - p)), NA)) # map to cluster
					s.to.c <- sapply(clusters$proportion, function(c) which.min(abs(cfi.s - c)))
					
					
					# Likelihood
					L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(altCount[hh],tumDepth[hh], ifelse(pp==1, 1-.Machine$double.eps, pp), rho=rho, xmin=pmin(altCount[hh],0)) + .Machine$double.eps), ncol=length(whichStates))
					
					# Power
					power.sm <- colMeans(matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) 1-ptrbetabinom(pmin(altCount[hh],xmin-1),tumDepth[hh],ifelse(pp==1, 1-.Machine$double.eps, pp), rho=rho, xmin=0)), ncol=length(whichStates)), na.rm=TRUE)

					if(globalIt==2){
						P.m.sX <- P[[h[i]]][,"P.m.sX"]
						power.s <- sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) # Power for state
						power.s[!is.na(c.to.s)] <- power.c[na.omit(c.to.s)]
						power.m.s <- power.sm # Relative power of m states (local) conditioned on s (global).
						for(s in unique(cnStates[whichStates,"s"])) power.m.s[cnStates[whichStates,"s"]==s] <- power.m.s[cnStates[whichStates,"s"]==s] / max(power.m.s[cnStates[whichStates,"s"]==s]) #power.s[s]
					}else{ # First iteration
						power.m.s <- power.sm
						power.s <- rep(1,length(whichStates))
					}
					
					mm <- function(x) {
						x <- factor(x)
						if(nlevels(x) > 1) t(model.matrix( ~ x + 0)) else matrix(1, ncol=length(x))
					}
					
					# EM algorithm (mixture fitting) for pi
					P.m.sX <- cnStates[whichStates,"pi.m.s"]
					s.from.m <- mm(cnStates[whichStates,"s"]) # indicator matrix to map
					for(em.it in 1:100){
						P.xsm <- L * rep(pi.s[cnStates[whichStates,"s"]] * P.m.sX / power.m.s / power.s[cnStates[whichStates,"s"]], each=nrow(L)) # P(X,s,m)
						P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
						P.sm.X <- pmax(.Machine$double.xmin,colMeans(P.sm.x)) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
						if(em.it==100) break
						P.s.X <- s.from.m %*% P.sm.X 
						P.m.sX <- P.sm.X / P.s.X[cnStates[whichStates,"s"]]
					}
					
					toTime <- function(cnStates, P.m.sX, s.m) {
						mAnc <- cnStates[,"m"] - cnStates[,"minDelta"] - cnStates[,"majDelta"]
						mAnc.s <- factor(paste(mAnc, cnStates[,"s"], sep="."))
						n <- (mAnc <= cnStates[,"majCNanc"]) + (mAnc <= cnStates[,"minCNanc"] )
						mAnc.s.from.m <- mm(x = mAnc.s)# indicator matrix to map
						return((mAnc.s.from.m[mAnc.s,] %*% P.m.sX) / (s.m[cnStates[,"s"],] %*% (P.m.sX * mAnc)) *  (cnStates[,"majCNanc"] + cnStates[,"minCNanc"]) / n)
					}
					
					T.m.sX <- toTime(cnStates, P.m.sX, s.from.m) 
					
					if(globalIt==1){
						p <- (sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) * nrow(L))[s.to.c]
						if(!any(is.na(p) | is.nan(p)))
							power.c <- power.c + p 
					}
					
					# Bootstrapping for CIs
					if(globalIt==2){
						b.m.sX <- if(n.boot>0) sapply(1:n.boot, function(foo){
												L <- rbind(L, rep(1e-3, each=ncol(L))) #add an uniformative row
												L <- L[sample(1:nrow(L), replace=TRUE),,drop=FALSE]
												P.m.sX <- cnStates[whichStates,"pi.m.s"]
												for(em.it in 1:50){
													P.xsm <- L * rep(pi.s[cnStates[whichStates,"s"]] * P.m.sX / power.m.s / power.s[cnStates[whichStates,"s"]], each=nrow(L)) # P(X,s,m)
													P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
													P.sm.X <- colMeans(P.sm.x) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
													P.s.X <- s.from.m %*% P.sm.X 
													P.m.sX <- P.sm.X / P.s.X[cnStates[whichStates,"s"]]
												}
												return(P.m.sX)
											}) else NA
						if(n.boot>0) try({
										CI.m.sX <- apply(b.m.sX, 1, quantile, c(0.025, 0.975))
										cnStates[,"P.m.sX.lo"] <- CI.m.sX[1,] 
										cnStates[,"P.m.sX.up"] <- CI.m.sX[2,]
										B.m.sX <- toTime(cnStates = cnStates, P.m.sX = b.m.sX, s.m = s.from.m)
										C.m.sX <- apply(B.m.sX, 1, quantile, c(0.025, 0.975))
										cnStates[,"T.m.sX.lo"] <- C.m.sX[1,] 
										cnStates[,"T.m.sX.up"] <- C.m.sX[2,]
									}, silent=TRUE)
					}
					
					P.sm.x[apply(is.na(P.sm.x)|is.nan(P.sm.x),1,any),] <- NA
					cnStates[,"P.m.sX"] <- P.m.sX
					cnStates[,"T.m.sX"] <- T.m.sX
					cnStates[,"power.s"] <- power.s[cnStates[whichStates,"s"]]
					cnStates[,"power.m.s"] <- power.m.s
					
				}else{
					cnStates <- bb$timing_param[[h[i]]]
					if(is.null(cnStates)) next
					P.sm.x <- posteriorMutCN(x=altCount[hh],n=tumDepth[hh], cnStates=cnStates, xmin=0, rho=rho)
				}
				
				# Tail probs
				pMutCNTail <- matrix(sapply(pmin(cnStates[,"f"],1), function(pp) ptrbetabinom(altCount[hh],tumDepth[hh], ifelse(pp==1, 1-.Machine$double.eps, pp), rho=rho, xmin=pmin(altCount[hh],xmin))), ncol=nrow(cnStates)) #%*% c(pi.s[cnStates[whichStates,"state"]] * P.m.sX)				
				
				P[[h[i]]] <- cnStates
				if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]
				
				w <- apply(P.sm.x, 1, function(x) if(any(is.na(x))) NA else which.max(x) )
				if(all(is.na(w))) next
				
				D[hh, "pSub"] <- rowSums(P.sm.x[, !cnStates[,"clonalFlag"], drop=FALSE])
				D[hh, "pGain"] <- rowSums(P.sm.x[, cnStates[,"clonalFlag"] & cnStates[,"m"] > 1.00001 + cnStates[,"majDelta"] + cnStates[,"minDelta"], drop=FALSE])
				#D[hh, "pSingle"] <- rowSums(P.sm.x[, cnStates[1:k,"state"] %in% which(clonalFlag) & cnStates[1:k,"m"]<=1, drop=FALSE])
				D[hh, "pSingle"] <-  1 - D[hh, "pSub"] - D[hh, "pGain"]			

				#currently doesn't work
				#D[hh,"pAllSubclones"] <- as(DataFrame(t(P.sm.x[, !cnStates[,"clonalFlag"], drop=FALSE])),"List")
				
				D[hh,"MutCN"]  <- cnStates[w,"m"]
				D[hh,"MutDeltaCN"]  <- cnStates[w,"majDelta"] + cnStates[w,"minDelta"]
				D[hh,"MinCN"] <- cnStates[w,"minCNanc"]
				D[hh,"MajCN"] <- cnStates[w,"majCNanc"]
				D[hh,"MinDerCN"] <- cnStates[w,"minCNder"]
				D[hh,"MajDerCN"] <- cnStates[w,"majCNder"]
				
				D[hh,"CNF"]  <- cnStates[w,"cfi"]
				D[hh,"pMutCN"] <- sapply(seq_along(w), function(i) P.sm.x[i,w[i]])
				D[hh,"pMutCNTail"] <- sapply(seq_along(w), function(i) pMutCNTail[i,w[i]])
				#D[hh,"altCount"] <- altCount[hh]
				#D[hh,"wtCount"] <- tumDepth[hh] - altCount[hh]
			}		
		}
		if(globalIt==1){
			power.c <- power.c / sum(!is.na(D[,"pSub"]))
			if(any(is.na(power.c) | power.c==0)) {
				break # Cancel 2nd iteration
				warning("Power calculation yielded NA, aborting.")
			}
			if(any(power.c > 1)) {
				warning("Calculated power > 1, rounding down.")
				power.c <- pmin(1, power.c)
			}
		}
	}
	return(list(D=D, P=P, power.c=power.c))
}

#' Compute mutation time from copy number gains and point mutations
#' @param vcf A vcf object of ssnms. See VariantAnnotation::readVcf()
#' @param cn The copy number as a GRanges() object, meta data in consensus format. See loadBB()
#' @param clusters A data.frame with the cluster proportion and n_ssms
#' @param purity The purity of the samples
#' @param gender 'male' or 'female'
#' @param isWgd TRUE/FALSE 
#' @param xmin min read support. Needed for power calculations
#' @param rho Dispersion parameter
#' @return A list with elements (V: Data.Frame with variant-specific timing information, can be added to vcf object; T: DataFrame with timing information to be added to CN Ranges; power.c power of each cluster).
#' @example inst/example/example.R
#' 
#' @author mg14
#' @export
mutationTime <- function(vcf, cn, clusters=data.frame(cluster=1, proportion=max(cn$clonal_frequency,na.rm=TRUE), n_ssms=100), purity=max(cn$clonal_frequency,na.rm=TRUE), gender='female', isWgd= FALSE, xmin=3, rho=0, n.boot=200){
	MT <- computeMutCn(vcf, cn, clusters, purity, gender, isWgd, xmin, rho, n.boot)
	MT$D$CLS <- classifyMutations(as.data.frame(MT$D))
	n.snv_mnv <- countOverlaps(cn,vcf)
	time <- mtToTime(cn, timing_param=MT$P, n.snv_mnv=n.snv_mnv)
	T <- DataFrame(timing_param=List(MT$P), time, n.snv_mnv=n.snv_mnv)
	return(list(V=MT$D, T=T))
}


#' Return header elements
#' @return DataFrame() to be added to VCF header
#' 
#' @author mg14
#' @export
mtHeader <- function() {
	DataFrame(Number=c(1,1,1,1,1,1,1,1,".",1,1,1,1,".",1),Type=c("Float","Float","Integer","Integer","Integer","Integer","Float","Float","Integer","Float","Float","Float","Float","String","String"), 
			Description=c("Mutation copy number","Change in MutCN between ancestral and derived state","Major copy number (ancestral)","Minor copy number (ancestral)","Major copy number (derived)","Minor copy number (derived)","Copy number frequency (relative to all cancer cells)", "MutCN probability","BB segment ids","Posterior prob: Early clonal","Posterior prob: Late clonal","Posterior prob: Subclonal", "Tail prob of mixture model", "Assignment probability of mutation to subclone","Mutation Time {clonal [early], clonal [late], clonal [NA], subclonal} - MAP assignments"),
			row.names=c("MutCN","MutDeltaCN","MajCN","MinCN","MajDerCN","MinDerCN","CNF","pMutCN","CNID","pGain","pSingle","pSub","pMutCNTail","pAllSubclones","CLS"))
}

mcnHeader <- function() {
	DataFrame(Number=c(1,1,1,1,1,1,1,1,".",1,1,1,1,"."),Type=c("Float","Float","Integer","Integer","Integer","Integer","Float","Float","Integer","Float","Float","Float","Float","String"), 
			Description=c("Mutation copy number","Change in MutCN between ancestral and derived state","Major copy number (ancestral)","Minor copy number (ancestral)","Major copy number (derived)","Minor copy number (derived)","Copy number frequency (relative to all cancer cells)", "MutCN probability","BB segment ids","Posterior prob: Early clonal","Posterior prob: Late clonal","Posterior prob: Subclonal", "Tail prob of mixture model", "Assignment probability of mutation to subclone"),
			row.names=c("MutCN","MutDeltaCN","MajCN","MinCN","MajDerCN","MinDerCN","CNF","pMutCN","CNID","pGain","pSingle","pSub","pMutCNTail","pAllSubclones"))
}

addMutCn <- function(vcf, bb=allBB[[meta(header(vcf))["ID",]]], clusters=allClusters[[meta(header(vcf))["ID",]]]){
  MCN = computeMutCn(vcf, bb, clusters)
	i = info(header(vcf))
	info(header(vcf)) <- rbind(i, mtHeader())
	info(vcf) <- cbind(info(vcf), MCN$D)
	return(vcf)
}

#' Add sample metadata entries to the VCF header
#' @param vcf A VCF object
#' @param purity Sample purity (Default: NULL)
#' @param ploidy Sample ploidy (Default: NULL)
#' @param sex Donor sex (Default: NULL)
#' @param is_wgd Sample whole genome doubling status (TRUE if a whole genome doubling has occurred) (Default: NULL)
#' @param subclonal_architecture Sample subclonal architecture data.frame (Default: NULL)
#' @param cluster_power Sample vector with a power estimate per mutation cluster (Default: NULL)
#' @return The VCF object with requested entries added to the metadata
#' @author sd11
addMetadataToHeader = function(vcf, purity=NULL, ploidy=NULL, sex=NULL, is_wgd=NULL, subclonal_architecture=NULL, cluster_power=NULL) {
  # convenience function that takes care of adding a metadata entry
  .addmeta = function(vcf, value, label) {
    temp = DataFrame(Value=value)
    rownames(temp) = label
    if ("META" %in% names(meta(header(vcf)))) {
      meta(header(vcf))[["META"]] = rbind(meta(header(vcf))[["META"]], temp)
    } else {
      meta(header(vcf))[["META"]] = temp
    }
    return(vcf)
  }
  
  if (!is.null(purity)) vcf = .addmeta(vcf, purity, "purity")
  if (!is.null(ploidy)) vcf = .addmeta(vcf, ploidy, "ploidy")
  if (!is.null(sex)) vcf = .addmeta(vcf, sex, "sex")
  if (!is.null(is_wgd)) vcf = .addmeta(vcf, is_wgd, "isWgd")
  if (!is.null(cluster_power)) vcf = .addmeta(vcf, paste(cluster_power, collapse=","), "clusterPower")
  
  # convenience function that collapses a data.frame column into a single string
  if (!is.null(subclonal_architecture)) {
    .collapsecol = function(subclonal_architecture, colname) paste(c(colname, ":", paste(subclonal_architecture[,colname], collapse=",")), collapse="")
    encoded_clusters = paste(sapply(colnames(subclonal_architecture), function(x) { .collapsecol(subclonal_architecture, x) }), collapse=";")
    vcf = .addmeta(vcf, encoded_clusters, "subclonalArchitecture")
  }
  
  return(vcf)
}

#' Convenience function to fetch info from metadata
#' @param vcf VCF object to query
#' @param label Row label to fetch from the metadata
#' @return The requested entry as a String, or NULL if the entry was not found
#' @author sd11
getMetaEntry = function(vcf, label) {
  if (label %in% rownames(meta(header(vcf))[["META"]])) {
    return(meta(header(vcf))[["META"]][label,"Value"])
  } else {
    return(NULL)
  }
}

#' Get sample purity from VCF metadata
#' @param VCF object to query
#' @return The sample purity, or NULL if the metadata entry is not available
#' @author sd11
#' @export
getPurity = function(vcf) {
  value = getMetaEntry(vcf, "purity")
  if (!is.null(value)) value = as.numeric(value)
  return(value)
}

#' Get sample ploidy from VCF metadata
#' @param VCF object to query
#' @return The sample ploidy, or NULL if the metadata entry is not available
#' @author sd11
#' @export
getPloidy = function(vcf) {
  value = getMetaEntry(vcf, "ploidy")
  if (!is.null(value)) value = as.numeric(value)
  return(value)
}

#' Get donor sex from VCF metadata
#' @param VCF object to query
#' @return The donor sex or, or NULL if the metadata entry is not available
#' @author sd11
#' @export
getSex = function(vcf) {
  value = getMetaEntry(vcf, "sex")
  if (!is.null(value)) value = as.character(value)
  return(value)
}

#' Get sample whole genome doubling status from VCF metadata
#' @param VCF object to query
#' @return The sample whole genome doubling status (TRUE if a WGD occurred, FALSE otherwise), or NULL if the metadata entry is not available
#' @author sd11
#' @export
getWGDstatus = function(vcf) {
  value = getMetaEntry(vcf, "isWgd")
  if (!is.null(value)) value = as.logical(value)
  return(value)
}

#' Get mutation cluster power from VCF metadata
#' @param VCF object to query
#' @return The mutation cluster power (vector with an entry per cluster), or NULL if the metadata entry is not available
#' @author sd11
#' @export
getClusterPower = function(vcf) {
  value = getMetaEntry(vcf, "clusterPower")
  if (!is.null(value)) value = as.numeric(unlist(strsplit(value, ",")))
  return(value)
}

#' Get sample subclonal architecture from VCF metadata
#' @param VCF object to query
#' @return A data.frame with the subclonal architecture of the sample, or NULL if the metadata entry is not available
#' @author sd11
#' @export
getSubclonalArchitecture = function(vcf) {
  value = getMetaEntry(vcf, "subclonalArchitecture")
  columns = lapply(unlist(strsplit(value, ";")), function(x) unlist(strsplit(x, ":")))
  columnnames = unlist(lapply(columns, function(x) x[1]))
  columnvalues = lapply(columns, function(x) unlist(strsplit(x[2], ",")))
  columnvalues = as.data.frame(do.call(cbind, columnvalues), stringsAsFactors=F)
  colnames(columnvalues) = columnnames
  
  for (columnname in c("cluster", "proportion", "n_ssms", "n_snvs", "ccf", "cp")) {
    if (columnname %in% colnames(columnvalues)) {
      columnvalues[,columnname] = as.numeric(columnvalues[,columnname])
    }
  }
  return(columnvalues)
}

#' Read a VCF file
#' 
#' This function overloads VariantAnnotation::readVcf as the pAllSubclones info column requires additional encoding to fit within the VCF standard
#' @param file Path to a VCF file
#' @return A VCF object
#' @importFrom VariantAnnotation readVcf
#' 
#' @author sd11
#' @export
readVcf = function(file, ...) {
  vcf = VariantAnnotation::readVcf(file, ...)
  if ("pAllSubclones" %in% colnames(info(vcf))) {
    info(vcf)$pAllSubclones = lapply(info(vcf)$pAllSubclones, function(x) as.numeric(unlist(strsplit(x, ","))))
  }
  return(vcf)
}

#' Write a VCF file
#' 
#' This function overloads VariantAnnotation::writeVcf as the pAllSubclones info column requires additional encoding to fit within the VCF standard
#' @param vcf A VCF object
#' @importFrom VariantAnnotation writeVcf
#' 
#' @author sd11
#' @export
writeVcf = function(vcf, ...) {
  if ("pAllSubclones" %in% colnames(info(vcf))) {
    info(vcf)$pAllSubclones = unlist(lapply(info(vcf)$pAllSubclones, function(x) paste(x, collapse=",")))
  }
  VariantAnnotation::writeVcf(vcf, ...)
}

#' Classify Mutations
#' @param x 
#' @param reclassify 
#' @return factor()
#' 
#' @importFrom S4Vectors as.matrix
#' @author mg14
classifyMutations <- function(x, reclassify=c("missing","all","none")) {
	reclassify <- match.arg(reclassify)
	if(nrow(x) ==0 )
		return(factor(NULL, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal")))
	if(class(x)=="CollapsedVCF")
		x <- info(x)
	.clsfy <- function(x) {
		cls <- x$CLS
		if(reclassify %in% c("missing", "none") &! is.null(cls)){
			if(all(unique(cls) %in% c("early", "late", "clonal", "subclonal")))
				cls <- factor(cls, levels=c("early", "late", "clonal", "subclonal"), labels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
			cls <- as.character(cls)
			cls[cls=="NA"] <- NA
			if(reclassify=="missing" & any(is.na(cls)))
				cls[is.na(cls)] <- paste(factor(apply(as.matrix(x[is.na(cls), c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
		}else{
			cls <- paste(factor(apply(as.matrix(x[, c("pGain","pSingle","pSub")]), 1, function(x) if(all(is.na(x))) NA else which.max(x)), levels=1:3, labels=c("clonal [early]", "clonal [late]","subclonal"))) ## reclassify missing
			
		}
		cls[x$pGain==0 & cls!="subclonal"] <- "clonal [NA]"
		if(!is.null(x$MajCN))
			cls[cls!="subclonal" & (x$MajCN == 1 | x$MinCN == 1) & abs(x$MutCN - x$MutDeltaCN -1) <= 0.0001] <- "clonal [NA]"
		cls <- factor(cls, levels=c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))
	}
	cls <- .clsfy(x)
	return(cls)
}

posteriorMutCN <- function(x,n, cnStates, xmin=3, rho=0.01){
	whichStates <- 1:nrow(cnStates)
	L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(x,n,pp, rho=rho, xmin=pmin(x,xmin)) + .Machine$double.eps), ncol=length(whichStates))
	P.xsm <- L * rep(cnStates[whichStates,"pi.s"] * cnStates[whichStates,"P.m.sX"] / cnStates[whichStates,"power.s"] / cnStates[whichStates,"power.m.s"], each=nrow(L)) # P(X,s,m)
	P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
	return(P.sm.x)
}

loadBB <- function(file){
  tab <- read.table(file, header=TRUE, sep='\t')
  GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])
}

pGainToTime <- function(vcf){
	P <- matrix(NA, nrow=nrow(vcf), ncol=4, dimnames=list(NULL, c("pEarly","pLate","pClonal[NA]","pSub")))
	P[,c("pEarly","pClonal[NA]","pSub")] <- as.matrix(info(vcf)[,c("pGain","pSingle","pSub")])
	biAllelicGain <- (info(vcf)$MajCN > 1 & (info(vcf)$MinCN > 1 | info(vcf)$MinCN == 0) & ! ((info(vcf)$MajCN == 1 | info(vcf)$MinCN == 1) & abs(info(vcf)$MutCN - info(vcf)$MutDeltaCN -1) <= 0.0001))
	w <- which(biAllelicGain)
	P[w, "pLate"] <- P[w, "pClonal[NA]"]
	P[w, "pClonal[NA]"] <- 0
	P[which(!biAllelicGain),"pLate"] <- 0
	return(P)
}

piToTime <- function(timing_param, type=c("Mono-allelic Gain","CN-LOH", "Bi-allelic Gain (WGD)")){
	type <- match.arg(type)
	evo <- NA
	w <- timing_param[,"s"]==1 &! timing_param[,"mixFlag"] ## Clonal states, not mixture (CN subclonal)
	M <- sum(w) ## Max multiplicity = Major copy number
	m <- sum(w & timing_param[,"n.m.s"]==2) ## Minor CN
	#evo <- paste0("1:1", paste0(2,":",m), if(M>2) paste0(3,":",) else "", collapse="->")
	t <- timing_param[(M:(M-1)),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-1
	if(nrow(t)==1) t <- rbind(t, NA)
	if(!any(is.na(t))){
		if(M==4) {
			if(timing_param[M-1,"P.m.sX.lo"] < 0.001){ ## No MCN 3
				if(type=="CN-LOH"){ ## Hotfix for 4+0, treated at 1:1 -> 2:0 -> 4:0
					t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
					evo <- "1:1->2:0->4:0"
				}
				else if(type=="Bi-allelic Gain (WGD)"){			
					if(m==2) {## Hotfix for 4+2 regions, treated at 1:1 -> 2:1 -> 4:2
						t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-2
						t[2,] <- pmax(0,2*t[2,] - t[1,])/3
						evo <- "1:1->2:1->4:2"
					} else if(m==4) {## Hotfix for 4+4 regions, treated at 1:1 -> 2:2 -> 4:4
						t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
						evo <- "1:1->2:2->4:4"
					}
				}			
			} else if(type=="Bi-allelic Gain (WGD)"){ ## Can't uniquely time second event
				t[2,] <- NA
			} 
			if(m==3) {
				t[2,] <- NA ## Don'time secondary 4+3 for now, needs more work
			} 
		} else {
			if(M==3 & type=="Bi-allelic Gain (WGD)") {## Hotfix for 3+2 regions, treated at 1:1 -> 2:2 -> 3:2
				t[2,] <- pmax(0,2*t[2,] - t[1,])
				evo <- "1:1->2:2->3:2"
			}
		}
	}
	colnames(t) <- c("","lo","up")
	t[,2] <- pmin(t[,1],t[,2])
	t[,3] <- pmax(t[,1],t[,3])
	t <- pmin(apply(t,2,cumsum),1) ## Times are actually deltas
	if(M < 3) t[2,] <- NA
	t[is.infinite(t)] <- NA
	rownames(t) <- c("",".2nd")
	return(c(t[1,],`.2nd`=t[2,])) ## Linearise
}

#' Convert timing parameters into timing estimates
#' @param cn Copy number input
#' @param timing_param 
#' @param pseudo.count 
#' @return data.frame()
#' 
#' @author mg14
#' @export
mtToTime <- function(cn, timing_param = cn$timing_param, n.snv_mnv = cn$n.snv_mnv, pseudo.count=5){
	sub <- duplicated(cn) 
	covrg <- countQueryHits(findOverlaps(cn, cn)) 
	maj <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "majCNanc"] else NA) #bb$major_cn
	min <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "minCNanc"] else NA) #bb$minor_cn
	type <- sapply(seq_along(cn), function(i){
				if(maj[i] < 2 | is.na(maj[i]) | sub[i] | (maj[i] > 4 & min[i] >= 2)) return(NA)
				type <- if(min[i]==1){ "Mono-allelic Gain" 
						}else if(min[i]==0){"CN-LOH"}
						else "Bi-allelic Gain (WGD)"
				return(type)
			})
	time <- t(sapply(seq_along(cn), function(i){
						if(sub[i] | is.na(type[i])) return(rep(NA,6)) 
						else piToTime(timing_param[[i]],type[i])
					}))
	
	res <- data.frame(type=factor(type, levels=c("Mono-allelic Gain","CN-LOH","Bi-allelic Gain (WGD)")), time=time)
	colnames(res) <- c("type","time","time.lo","time.up","time.2nd","time.2nd.lo","time.2nd.up")
	
	# posthoc adjustment of CI's
	res$time.up <- (pseudo.count + n.snv_mnv * res$time.up)/(pseudo.count + n.snv_mnv)
	res$time.lo <- (0 + n.snv_mnv * res$time.lo)/(pseudo.count + n.snv_mnv)
	res$time.2nd.up <- (pseudo.count + n.snv_mnv * res$time.2nd.up)/(pseudo.count + n.snv_mnv)
	res$time.2nd.lo <- (0 + n.snv_mnv * res$time.2nd.lo)/(pseudo.count + n.snv_mnv)
	
	res$time.star <- factor((covrg == 1) + (min < 2 & maj <= 2 | min==2 & maj==2) * (covrg == 1), levels=0:2, labels = c("*","**","***")) ## ***: 2+0, 2+1, 2+2; **: 2<n+1, {3,4}+2; *: subclonal gains
	res$time.star[is.na(res$time)] <- NA
	return(res)
}

averagePloidy <- function(bb) {
	c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
	sum(width(bb) * c * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

#' @importFrom VGAM rbetabinom
simulateMutations <- function(cn, purity=max(cn$clonal_frequency, na.rm=TRUE),  n=40, rho=0.01, xmin=3){
	g <- (averagePloidy(cn)*purity + 2*(1-purity))
	V <- list(VRanges())#VRanges()
	for(i in which(!duplicated(cn)))
		if(cn$n.snv_mnv[i]>1 & !is.null( cn$timing_param[[i]]))try({
						cnStates <- cn$timing_param[[i]]
						p <- cnStates[,"pi.s"]* if(!any(is.na(cnStates[,"P.m.sX"]))) cnStates[,"P.m.sX"] else cnStates[,"pi.m.s"]
						pwr <- cnStates[,"power.m.s"]#(cnStates[,"power.s"] * cnStates[,"power.m.s"])
						s <- sample(1:nrow(cnStates), size=pmax(1,ceiling(cn$n.snv_mnv[i] * (p %*% (1/pwr)))), prob=p, replace=TRUE)
						f <- cnStates[s,"f"]
						mu.c <- ((cn$major_cn[i] + cn$minor_cn[i])*purity + 2*(1-purity))/g * n
						c <- rnbinom(length(f), size=1/rho, mu=mu.c)
						x <- rbetabinom(n=length(f), size=c, prob=f, rho=rho)
						pos <- round(runif(length(f), min=start(cn)[i], max=end(cn)[i]))
						w <- which(x>=xmin)
						V[[i]] <- VRanges(seqnames=seqnames(cn)[i], IRanges(pos, width=1), ref="N", alt="A", totalDepth=c, altDepth=x)[w]
					})
	V <- do.call("c", V[!sapply(V, is.null)])
	sampleNames(V) <- "SAMPLE"
	v <- as(V, "VCF")
	info(header(v)) <- rbind(info(header(v)),DataFrame(Number=1, Type=rep("Integer",2),Description=c("Tumour ref count","Tumour alt count"), row.names=c("t_ref_count","t_alt_count")))
	info(v)$t_alt_count <- altDepth(V)
	info(v)$t_ref_count <- totalDepth(V) - altDepth(V)
	return(v)
}


.plotBB <- function(bb, ylim=c(0,max(max(bb$total_cn, na.rm=TRUE))), col=RColorBrewer::brewer.pal(4,"Set2"), type=c("lines","bars"), legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, xlim=c(min(chrOffset[as.character(seqnames(bb))]+start(bb)),max(chrOffset[as.character(seqnames(bb))]+end(bb)))){
	type <- match.arg(type)
	s <- c(1:22, "X","Y")
	l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
	names(l) <- s
	plot(NA,NA, ylab="Copy number",xlab="",xlim=xlim, ylim=ylim, xaxt="n")
	c <- cumsum(l)
	axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
	if(xaxt) mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
	#abline(v=c, lty=3)
	if(type=="lines"){
		x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
		x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
		lwd <- 5* bb$clonal_frequency / max(bb$clonal_frequency)
		segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=col[1], lwd=lwd)
		segments(x0=x0, bb$minor_cn -.125,x1, bb$minor_cn-.125, col=col[2], lwd=lwd)
		segments(x0=x0, bb$total_cn+.125,x1, bb$total_cn+.125, col=1, lwd=lwd)
#	cv <- coverage(bb)
#	cv <- cv[s[s%in%names(cv)]]
#	par(xpd=NA)
#	for(n in names(cv)){
#		cc <- cv[[n]]
#		segments(start(cc) + cumsum(l)[n] - l[n] ,-runValue(cc)/2,end(cc)+ cumsum(l)[n] - l[n], -runValue(cc)/2, col=4)
#	}
	}else{
		ub <- unique(bb)
		f <- findOverlaps(ub,bb)
		m <- t(model.matrix( ~ 0 + factor(queryHits(f))))
		ub$major_cn <- m %*% mg14::na.zero(bb$major_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
		ub$minor_cn <- m %*% mg14::na.zero(bb$minor_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
		ub$total_cn <- ub$major_cn + ub$minor_cn
		ub$clonal_frequency <- max(bb$clonal_frequency)
		x0 <- start(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
		x1 <- end(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
		rect(x0,0,x1, ub$major_cn, col=col[2], lwd=NA)
		rect(x0,ub$major_cn,x1, ub$total_cn, col=col[1], lwd=NA)
		abline(h = 1:floor(ylim[2]), lty=lty.grid, col=col.grid)
	}
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	if(legend){
		if(type=="lines") legend("topleft", c("Total CN","Major CN","Minor CN"), col=c("black", col[1:2]), lty=1, lwd=2, bg='white')
		else legend("topleft", c("Major CN","Minor CN"), fill=col[1:2], bg='white')
	}
}

.plotVcf <- function(vcf, bb, col = RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)], ID = meta(header(vcf))[[1]]["ID",1], legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, pch=16, pch.out=pch, cex=0.66, xlim=c(0,chrOffset["MT"])) {
	cls <- factor(paste(as.character(info(vcf)$CLS)), levels = c("clonal [early]","clonal [late]","clonal [NA]","subclonal" , "NA"))
	plot(NA,NA, xlab='', ylab="VAF", ylim=c(0,1), xlim=xlim, xaxt="n", cex=cex)
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	for(i in which(!sapply(bb$timing_param, is.null))) {
		s <- start(bb)[i]
		e <- end(bb)[i]
		x <- chrOffset[as.character(seqnames(bb)[i])]
		y <- bb$timing_param[[i]][,"f"]
		l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
		l[is.na(l)] <- 0
		if(any(is.na(c(s,e,x,y,l)))) next
		segments(s+x,y,e+x,y, lwd=l*4+.1)
		#text(x=(s+e)/2 +x, y=y, paste(signif(bb$timing_param[[i]][,"m"],2),signif(bb$timing_param[[i]][,"cfi"]/purityPloidy[meta(header(vv))["ID",1],"purity"],2), sep=":"), pos=3, cex=0.5)
	}
	points(start(vcf) + chrOffset[as.character(seqnames(vcf))], getAltCount(vcf)/getTumorDepth(vcf),col=col[cls],  pch=ifelse(info(vcf)$pMutCNTail < 0.025 | info(vcf)$pMutCNTail > 0.975, pch.out , pch),  cex=cex)				
	if(legend) legend("topleft", pch=19, col=col, legend=paste(as.numeric(table(cls)), levels(cls)), bg='white')
}

.timeToBeta <- function(time){
	mu <- time[,1]
	#if(any(is.na(time))) return(c(NA,NA))
	mu <- pmax(1e-3, pmin(1 - 1e-3, mu))
	v <- (0.5 * (pmax(mu,time[,3])-pmin(mu,time[,2])))^2
	alpha <- mu * (mu * (1-mu) / v - 1)
	beta <- (1-mu) *  (mu * (1-mu) / v - 1)
	return(cbind(alpha, beta))
}

.plotTiming <- function(bb, time=mcols(bb), col=paste0(RColorBrewer::brewer.pal(5,"Set2")[c(3:5)],"88"), legend=TRUE, col.grid='grey', lty.grid=1, xlim=c(0,chrOffset["MT"]), plot=2){
	plot(NA,NA, xlab='', ylab="Time [mutations]", ylim=c(0,1), xlim=xlim, xaxt="n")
	if(any(!is.na(bb$time)))
		tryCatch({
					bb <- bb[!is.na(bb$time)]
					s <- start(bb)
					e <- end(bb)
					x <- chrOffset[as.character(seqnames(bb))]
					y <- time[,"time"]
					rect(s+x,time[,"time.lo"],e+x,time[,"time.up"], border=NA, col=col[time[,"type"]], angle = ifelse(bb$time.star=="*" | is.na(bb$time.star),45,135), density=ifelse(bb$time.star == "***", -1, 72))
					segments(s+x,y,e+x,y)
					
					if("time.2nd" %in% colnames(time)){ 
						w <- !is.na(time[,"time.2nd"])
						if(sum(w) != 0 & plot==2){
							s <- start(bb)[w]
							e <- end(bb)[w]
							x <- chrOffset[as.character(seqnames(bb))][w]
							y <- time[w,"time.2nd"]
							rect(s+x,time[w,"time.2nd.lo"],e+x,time[w,"time.2nd.up"], border=NA, col=sub("88$","44",col)[as.numeric(time[w,"type"])], angle = ifelse(bb$time.star[w]=="*" | is.na(bb$time.star[w]),45,135), density=ifelse(bb$time.star[w] == "***", -1, 72))
							segments(s+x,y,e+x,y)
						}
					}
				}, error=function(x) warning(x))
	abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
	s <- c(1:22, "X","Y")
	l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
	names(l) <- s
	c <- cumsum(l)
	axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
	mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
	if(legend) legend("topleft", levels(time[,"type"]), fill=col, bg="white")
}

.betaFromCi <- function(x, weight.mode=5){
	if(any(is.na(x))) return(rep(NA,2))
	f <- function(par,x) {
		beta <- exp(par)
		sum((qbeta(c(0.025,0.975), beta[1], beta[2])-x[-1])^2)+weight.mode*((beta[1]-1)/(beta[1]+beta[2]-2)-x[1])^2
	}
	tryCatch(exp(optim(c(0.1,0.1), fn=f,x=x)$par), error=c(1,1))
}

.histBeta <- function(bb, time="time",n.min=10, s=seq(0.005,0.995,0.01)){
	s <- pmax(0.001,pmin(0.999, s))
	cols <- paste0(time,c("",".lo",".up"))
	w <- which(bb$n.snv_mnv > n.min & !is.na(mcols(bb)[cols[1]]) & !duplicated(bb))
	if(length(w)==0) return(rep(NA, length(s)))
	d <- apply(as.matrix(mcols(bb)[w,c(cols, "n.snv_mnv")]), 1, function(x){
				beta <- .betaFromCi(x[1:3])
				beta <- (beta * x[4] + 5*c(1,1))/(x[4]+5) # "posterior" with prior B(1,1)
				dbeta(s,beta[1],beta[2])
			})
	wd <- as.numeric(width(bb)[w])
	o <- d %*% wd
}



#' Plot timing
#' @param vcf The `VCF` file to plot, with timing information added 
#' @param cn 
#' @param sv 
#' @param title 
#' @param regions 
#' @param ylim.bb 
#' @param layout.height 
#' @param y1 
#' @return NULL
#' 
#' @author mg14
#' @export
plotSample <- function(vcf, cn, sv=NULL, title="", regions=NULL, ylim.bb=c(0,5), layout.height=c(4,1.2,3.5), y1=ylim.bb[2]-1) {
	if(is.null(regions)) regions <- refLengths[1:24]
	p <- par()
	layout(matrix(1:3, ncol=1), height=layout.height)
	par(mar=c(0.5,3,0.5,0.5), mgp=c(2,0.25,0), bty="L", las=2, tcl=-0.25, cex=1)
	xlim=c(min(chrOffset[as.character(seqnames(regions))]+start(regions)),max(chrOffset[as.character(seqnames(regions))]+end(regions)))
	bbb <- cn[cn %over% regions]
	.plotVcf(vcf[vcf %over% regions], bbb, legend=FALSE, col.grid='white',  xaxt=FALSE, cex=0.33, xlim=xlim)
	mtext(line=-1, side=3, title, las=1)
	.plotBB(bbb, ylim=ylim.bb, legend=FALSE, type='bar', col.grid='white', col=c("lightgrey", "darkgrey"), xaxt=FALSE, xlim=xlim)
	tryCatch({
				par(xpd=NA)
				if(!is.null(sv))
					.plotSv(sv, y1=y1, regions=regions, add=TRUE)
				par(xpd=FALSE)
			}, error=function(x) warning(x))
	par(mar=c(3,3,0.5,0.5))
	.plotTiming(bbb, xlim=xlim, legend=FALSE, col.grid=NA)
	if(length(regions) == 1)
		axis(side=1, at=pretty(c(start(regions), end(regions)))+chrOffset[as.character(seqnames(regions))], labels=sitools::f2si(pretty(c(start(regions), end(regions)))))
	if(any(!is.na(cn$time))){
		y0 <- seq(0.005,0.995,0.01)
		s <- .histBeta(cn)
		g <- colorRampPalette(RColorBrewer::brewer.pal(4,"Set1")[c(3,2,4)])(100)
		segments(x0=chrOffset["MT"] ,y0=y0,x1=chrOffset["MT"] + s/max(s) * 1e8, col=g, lend=3)
		getMode <- function(s){
			if(all(is.na(s))) return(NA)
			w <- which.max(s)
			if(w %in% c(1, length(s))){
				m <- which(c(0,diff(s))>0 & c(diff(s),0)<0)
				if(length(m)==0) return(w)
				m <- m[which.max(s[m])]
				return(if(s[w] > 2*s[m]) w else m) 
			} else return(w)
		}
		abline(h=y0[getMode(s)], lty=5)
		if("time.2nd" %in% colnames(mcols(cn))) if(any(!is.na(cn$time.2nd))){
				s2 <- .histBeta(cn, time="time.2nd")
				segments(x0=chrOffset["MT"] + s/max(s) * 1e8 ,y0=y0,x1=chrOffset["MT"] + s/max(s) * 1e8 + s2/max(s) * 1e8, col=paste0(g,"44"), lend=3)
				abline(h=y0[getMode(s2)], lty=3)
				
			}
	}
	#print(w)
	par(p[setdiff(names(p), c("cin","cra","csi","cxy","din","page"))])
}

.plotSv <- function(sv, y0=0,y1=y0, h=1, col=paste0(RColorBrewer::brewer.pal(5,"Set1"),"44"), regions=refLengths[1:24], add=FALSE){
	if(add==FALSE){
		s <- c(1:22, "X","Y")
		l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
		names(l) <- s
		plot(NA,NA, ylab="Copy number",xlab="",xlim=xlim, ylim=ylim, xaxt="n")
		c <- cumsum(l)
		axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
		if(length(regions) == 1)
			axis(side=1, at=pretty(c(start(regions), end(regions)))+chrOffset[as.character(seqnames(regions))], labels=sitools::f2si(pretty(c(start(regions), end(regions)))))
		if(xaxt) mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
	}
	#r <- rowRanges(sv)
	#a <- unlist(alt(sv))
	vs <- GRanges(info(sv)$MATECHROM, IRanges(info(sv)$MATEPOS, width=1))
	l <- 20
	x0 <- seq(0,1,l=l) 
	y2 <- x0*(1-x0)*4*h
	cls <- factor(as.character(info(sv)$SVCLASS), levels=c("DEL", "DUP", "h2hINV","t2tINV","TRA"))
	w <- which(sv %over% regions | vs %over% regions)
	for(i in w)
		try({
					x <- seq(start(sv)[i] + chrOffset[as.character(seqnames(sv)[i])], start(vs)[i] + chrOffset[as.character(seqnames(vs)[i])], length.out=l)
					x <- c(x[1], x, x[length(x)])
					y <- y1 + y2 * if(grepl("INV", cls[i])) -1 else 1
					y <- c(y0, y , y0)
					lines(x, y, col=col[cls[i]])
					#segments(x0=c(x[1], x[l]), x1=c(x[1],x[l]), y0=y0, y1=y1, col=col[cls[i]])
				})
}

