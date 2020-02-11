cn <- MutationTimeR:::refLengths[1:23]
cn$major_cn <- 3
cn$minor_cn <- 1
cn$clonal_frequency <- 0.8
cn$timing_param <- MutationTimeR:::defineMcnStates(cn, purity=0.8, clusters=data.frame(cluster=1, proportion=0.8, n_ssms=100))
for(i in seq_along(cn)){
	cn$timing_param[[i]][,"P.m.sX"] <- c(0.6,0.1,0.3)
	cn$timing_param[[i]][,"power.m.s"] <- rep(1,3)
}
cn$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6

vcf <- MutationTimeR:::simulateMutations(cn)
cn_timed <- cn[,c("major_cn","minor_cn","clonal_frequency")]

mt <- mutationTime(vcf, cn_timed, n.boot=10)
mcols(cn_timed) <- cbind(mcols(cn_timed),mt$T)

info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)

plotSample(vcf,cn_timed)
