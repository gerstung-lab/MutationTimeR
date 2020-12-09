set.seed(42)
cn <- MutationTimeR:::refLengths[1:23]
t1 <- 0.7
t2 <- 0.9
purity <- 0.6
clusters <- data.frame(cluster=1:3, proportion=c(purity,purity * 0.8, purity * 0.4), n_ssms=c(100,20, 50))
time2pi <- function(N,n,t1,t2){
	if(N==2 & n ==1)
		pi <-  c(3-2*t1, t1)
	else if(N==2 & n %in% c(0,2))
		pi <- c(2 -2*t1, t1)
	else pi <- 1
	pi <- pi/sum(pi)
}
cn$major_cn <- 2 #sample(1:2, 23, replace=TRUE)
cn$minor_cn <- sample(0:2, 23, replace=TRUE)
cn$clonal_frequency <- purity
cn$timing_param <- MutationTimeR:::defineMcnStates(cn, purity=purity, clusters=clusters, deltaFreq = 0)
for(i in seq_along(cn)){
	pi <- time2pi(cn$major_cn[i], cn$minor_cn[i], t1, t2)
	pi_sub <- clusters$n_ssms
	pi_sub[1] <- pi_sub[1] * (t1*2 + (1-t1)*(cn$major_cn[i]+ cn$minor_cn[i])) / 2
	pi_sub[-1] <- pi_sub[-1] * (cn$major_cn[i]+ cn$minor_cn[i]) / 2
	pi_sub <- pi_sub/sum(pi_sub)
	#print(cbind(pi_sub, clusters$n_ssms/sum(clusters$n_ssms)))
	pi <- c(pi * pi_sub[1], pi_sub[-1])
	cn$timing_param[[i]][,"P.m.sX"] <- pi
	cn$timing_param[[i]][,"power.m.s"] <- rep(1, length(pi))
}
cn$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6 * 10

vcf <- MutationTimeR:::simulateMutations(cn, rho=0.01, n=40)
cn_timed <- cn[,c("major_cn","minor_cn","clonal_frequency")]

clusters_o <- clusters#[-2,]
bias <- 0.9
clusters_o$proportion <- clusters_o$proportion * bias
cn_timed$clonal_frequency <- cn_timed$clonal_frequency * bias
#clusters_o$proportion[2] <- 0.2
#clusters_o$n_ssms[2] <- 20

mt <- mutationTime(vcf, cn_timed, clusters=clusters_o, n.boot=10)
mcols(cn_timed) <- cbind(mcols(cn_timed),mt$T)

info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)

plotSample(vcf,cn_timed)
