#' # Calculation of mutation copy numbers and timing parameters
#' Minimal example
#' source("/Users/mg14/Projects/PCAWG-11/code/MutationTime.R")
#' vcf <- readVcf("final/final_consensus_12oct_passonly/snv_mnv/0040b1b6-b07a-4b6e-90ef-133523eaf412.consensus.20160830.somatic.snv_mnv.vcf.gz",genome="GRCh37")
#' bb <- loadBB("dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/4_copynumber/0040b1b6-b07a-4b6e-90ef-133523eaf412_segments.txt.gz")
#' clusters = read.table("dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/2_subclones/0040b1b6-b07a-4b6e-90ef-133523eaf412_subclonal_structure.txt.gz", header=TRUE, sep="\t")
#' purityPloidy <- read.table("dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/1_purity_ploidy/purity_ploidy.txt, header=TRUE, sep="\t")
#' purity <- purityPloidy[2,2]
#' MCN <- computeMutCn(vcf[1:1000], bb, clusters, purity)
#' # Check output
#' # Mutation Annotation
#' head(MCN$D)
#' # Classify as basic clonal states
#' table(classifyMutations(MCN$D))
#' # Timing parameters
#' MCN$P[[1]]
#' # Extract timing of segments
#' bb$timing_param <- MCN$P
#' bbToTime(bb)
#' PRE: Purity must be equal to the main clone, otherwise spurious results can be generated
#' PRE: CNV segments don't have to overlap, otherwise an error is generated: "Error in if (!mixFlag[j]) { : missing value where TRUE/FALSE needed"

require(VariantAnnotation)
require(VGAM)



#' Extract tumour depth
#' @param vcf 
#' @return 
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
#' @return 
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

dtrbinom <- function(x, size, prob, xmin=0) dbinom(x,size,prob) / pbinom(xmin-1, size, prob, lower.tail=FALSE)
pbetabinom <- function(x, size, prob, rho){
	if(rho!=0)
		VGAM::pbetabinom(x, size, prob, rho)
	else
		pbinom(x, size, prob)
}
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
#' @param bb Copy number with consensus annotation meta data
#' @param clusters 
#' @param purity 
#' @param gender 
#' @param isWgd 
#' @return list of length nrow(bb), can be added to mcols(bb)
#' 
#' @author mg14
#' @export
defineMcnStates <- function(bb, clusters, purity, gender='female', isWgd= FALSE){
	P <- vector(mode='list', length(bb))
	uniqueBB <- unique(bb)
	overlaps <- findOverlaps(uniqueBB, bb)
	
	majorCN <- split(bb$major_cn[subjectHits(overlaps)], queryHits(overlaps))
	m <- bb$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
	if(gender=='male')
		m[as.character(seqnames(bb)) %in% c('X','Y')] <- 0
	minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps))	
	h <- selectHits(overlaps, "first")
	H <- selectHits(overlaps, "last")
	
	cnNormal <- 2 - (gender=='male' & seqnames(uniqueBB)=="X" | seqnames(uniqueBB)=="Y")
	
	# Fix cluster and purity discrepancies
	clusters$proportion[which.max(clusters$proportion)] <- purity
	
	cloneFreq <- split(bb$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
	cnStates <- matrix(0, nrow=10000, ncol=6)
	colnames(cnStates) <- c("state","m","f","n.m.s","pi.m.s","s")
	
	power.c <- rep(0, nrow(clusters))
	
	deltaFreq <- 0.05 # merge clusters withing deltaFreq
	
	
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


#' Compute Mutation copy numbers
#' @param vcf A vcf object of ssnms. See VariantAnnotation::readVcf()
#' @param bb The copy number as a GRanges() object, meta data in consensus format. See loadBB()
#' @param clusters A data.frame with the cluster proportion and n_ssms
#' @param purity The purity of the samples
#' @param gender 'male' or 'female'
#' @param isWgd TRUE/FALSE 
#' @param xmin min read support. Needed for power calculations
#' @param rho Dispersion parameter
#' @return A list with elements (D: Annotated Data.Frame, can be added to vcf object; P: Timing parameters to be added to CN Ranges; power.c power of each cluster).
#' 
#' @author mg14
#' @export
computeMutCn <- function(vcf, bb, clusters, purity, gender='female', isWgd= FALSE, xmin=3, rho=0, n.boot=200){
	n <- nrow(vcf)
	D <- DataFrame(MutCN=rep(NA,n), MutDeltaCN=rep(NA,n), MajCN=rep(NA,n), MinCN=rep(NA,n), MajDerCN=rep(NA,n), MinDerCN=rep(NA,n), CNF=rep(NA,n), CNID =as(vector("list", n),"List"), pMutCN=rep(NA,n), pGain=rep(NA,n),pSingle=rep(NA,n),pSub=rep(NA,n), pMutCNTail=rep(NA,n))	
	P <- defineMcnStates(bb,clusters, purity, gender, isWgd)
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
					L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(altCount[hh],tumDepth[hh],pp, rho=rho, xmin=pmin(altCount[hh],0)) + .Machine$double.eps), ncol=length(whichStates))
					
					# Power
					power.sm <- colMeans(matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) 1-ptrbetabinom(pmin(altCount[hh],xmin-1),tumDepth[hh],pp, rho=rho, xmin=0)), ncol=length(whichStates)), na.rm=TRUE)
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
				pMutCNTail <- matrix(sapply(pmin(cnStates[,"f"],1), function(pp) ptrbetabinom(altCount[hh],tumDepth[hh],pp, rho=rho, xmin=pmin(altCount[hh],xmin))), ncol=nrow(cnStates)) #%*% c(pi.s[cnStates[whichStates,"state"]] * P.m.sX)				
				
				P[[h[i]]] <- cnStates
				if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]
				
				w <- apply(P.sm.x, 1, function(x) if(any(is.na(x))) NA else which.max(x) )
				if(all(is.na(w))) next
				
				D[hh, "pSub"] <- rowSums(P.sm.x[, !cnStates[,"clonalFlag"], drop=FALSE])
				D[hh, "pGain"] <- rowSums(P.sm.x[, cnStates[,"clonalFlag"] & cnStates[,"m"] > 1.00001 + cnStates[,"majDelta"] + cnStates[,"minDelta"], drop=FALSE])
				#D[hh, "pSingle"] <- rowSums(P.sm.x[, cnStates[1:k,"state"] %in% which(clonalFlag) & cnStates[1:k,"m"]<=1, drop=FALSE])
				D[hh, "pSingle"] <-  1 - D[hh, "pSub"] - D[hh, "pGain"]			
				
				D[hh,"MutCN"]  <- cnStates[w,"m"]
				D[hh,"MutDeltaCN"]  <- cnStates[w,"majDelta"] + cnStates[w,"minDelta"]
				D[hh,"MinCN"] <- cnStates[w,"minCNanc"]
				D[hh,"MajCN"] <- cnStates[w,"majCNanc"]
				D[hh,"MinDerCN"] <- cnStates[w,"minCNder"]
				D[hh,"MajDerCN"] <- cnStates[w,"majCNder"]
				
				D[hh,"CNF"]  <- cnStates[w,"cfi"]
				D[hh,"pMutCN"] <- sapply(seq_along(w), function(i) P.sm.x[i,w[i]])
				D[hh,"pMutCNTail"] <- sapply(seq_along(w), function(i) pMutCNTail[i,w[i]])
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
	return(list(D=D,P=P, power.c=power.c))
}

mcnHeader <- function() {
	DataFrame(Number=c(1,1,1,1,1,1,1,1,".",1,1,1,1),Type=c("Float","Float","Integer","Integer","Integer","Integer","Float","Float","Integer","Float","Float","Float","Float"), 
			Description=c("Mutation copy number","Change in MutCN between ancestral and derived state","Major copy number (ancestral)","Minor copy number (ancestral)","Major copy number (derived)","Minor copy number (derived)","Copy number frequency (relative to all cancer cells)", "MutCN probability","BB segment ids","Posterior prob: Early clonal","Posterior prob: Late clonal","Posterior prob: Subclonal", "Tail prob of mixture model"),
			row.names=c("MutCN","MutDeltaCN","MajCN","MinCN","MajDerCN","MinDerCN","CNF","pMutCN","CNID","pGain","pSingle","pSub", "pMutCNTail"))
}

addMutCn <- function(vcf, bb=allBB[[meta(header(vcf))["ID",]]], clusters=allClusters[[meta(header(vcf))["ID",]]]){
	i = info(header(vcf))
	info(header(vcf)) <- rbind(i, mcnHeader())
	info(vcf) <- cbind(info(vcf), computeMutCn(vcf, bb, clusters)$D)
	return(vcf)
}

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
	w <- timing_param[,"s"]==1 &! timing_param[,"mixFlag"] 
	n <- sum(w)
	t <- timing_param[n,c("T.m.sX","T.m.sX.lo","T.m.sX.up")]
	names(t) <- c("","lo","up")
	t[2] <- min(t[1],t[2])
	t[3] <- max(t[1],t[3])
	return(pmin(t,1))
}

bbToTime <- function(bb, timing_param = bb$timing_param, pseudo.count=5){
	sub <- duplicated(bb) 
	covrg <- countQueryHits(findOverlaps(bb, bb)) 
	maj <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "majCNanc"] else NA) #bb$major_cn
	min <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "minCNanc"] else NA) #bb$minor_cn
	type <- sapply(seq_along(bb), function(i){
				if(maj[i] < 2 | is.na(maj[i]) | sub[i] | (maj[i] > 2 & min[i] >= 2)) return(NA)
				type <- if(min[i]==1){ "Mono-allelic Gain" 
						}else if(min[i]==0){"CN-LOH"}
						else "Bi-allelic Gain (WGD)"
				return(type)
			})
	time <- t(sapply(seq_along(bb), function(i){
						if(sub[i] | is.na(type[i])) return(c(NA,NA,NA)) 
						else piToTime(timing_param[[i]],type[i])
					}))
	
	res <- data.frame(type=factor(type, levels=c("Mono-allelic Gain","CN-LOH","Bi-allelic Gain (WGD)")), time=time)
	colnames(res) <- c("type","time","time.lo","time.up")
	
	# posthoc adjustment of CI's
	res$time.up <- (pseudo.count + bb$n.snv_mnv * res$time.up)/(pseudo.count + bb$n.snv_mnv)
	res$time.lo <- (0 + bb$n.snv_mnv * res$time.lo)/(pseudo.count + bb$n.snv_mnv)
	res$time.star <- factor((covrg == 1) + (min < 2 & maj <= 2 | min==2 & maj==2) * (covrg == 1), levels=0:2, labels = c("*","**","***"))
	res$time.star[is.na(res$time)] <- NA
	return(res)
}

simulateMutations <- function(bb, purity=max(bb$clonal_frequency, na.rm=TRUE),  n=40, rho=0.01, xmin=3){
	g <- (averagePloidy(bb)*purity + 2*(1-purity))
	V <- list(VRanges())#VRanges()
	for(i in which(!duplicated(bb)))
		if(bb$n.snv_mnv[i]>1 & !is.null( bb$timing_param[[i]]))try({
						cnStates <- bb$timing_param[[i]]
						p <- cnStates[,"pi.s"]* if(!any(is.na(cnStates[,"P.m.sX"]))) cnStates[,"P.m.sX"] else cnStates[,"pi.m.s"]
						pwr <- cnStates[,"power.m.s"]#(cnStates[,"power.s"] * cnStates[,"power.m.s"])
						s <- sample(1:nrow(cnStates), size=pmax(1,ceiling(bb$n.snv_mnv[i] * (p %*% (1/pwr)))), prob=p, replace=TRUE)
						f <- cnStates[s,"f"]
						mu.c <- (bb$total_cn[i]*purity + 2*(1-purity))/g * n
						c <- rnbinom(length(f), size=1/rho, mu=mu.c)
						x <- rbetabinom(n=length(f), size=c, prob=f, rho=rho)
						pos <- round(runif(length(f), min=start(bb)[i], max=end(bb)[i]))
						w <- which(x>=xmin)
						V[[i]] <- VRanges(seqnames=seqnames(bb)[i], IRanges(pos, width=1), ref="N", alt="A", totalDepth=c, altDepth=x)[w]
					})
	V <- do.call("c", V[!sapply(V, is.null)])
	sampleNames(V) <- "SAMPLE"
	v <- as(V, "VCF")
	info(header(v)) <- rbind(info(header(v)),info(header(finalSnv[[1]]))[c("t_ref_count","t_alt_count"),])
	info(v)$t_alt_count <- altDepth(V)
	info(v)$t_ref_count <- totalDepth(V) - altDepth(V)
	return(v)
}

