set.seed(42)
cn <- MutationTimeR:::refLengths[1:23]
t1 <- 0.7
t2 <- 0.9
purity <- 0.7
clusters <- data.frame(cluster=1:2, proportion=c(purity,purity * 0.4), n_ssms=c(100,0.01))
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
cn$timing_param <- MutationTimeR:::defineMcnStates(cn, purity=purity, clusters=clusters)
for(i in seq_along(cn)){
	pi <- time2pi(cn$major_cn[i], cn$minor_cn[i], t1, t2)
	pi_sub <- clusters$n_ssms
	pi_sub <- pi_sub/sum(pi_sub)
	pi <- c(pi * pi_sub[1], pi_sub[2])
	cn$timing_param[[i]][,"P.m.sX"] <- pi
	cn$timing_param[[i]][,"power.m.s"] <- rep(1, length(pi))
}
cn$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6 * 10

vcf <- MutationTimeR:::simulateMutations(cn)
cn_timed <- cn[,c("major_cn","minor_cn","clonal_frequency")]

mt <- mutationTime(vcf, cn_timed, clusters=clusters, n.boot=10)
mcols(cn_timed) <- cbind(mcols(cn_timed),mt$T)

info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)

plotSample(vcf,cn_timed)
