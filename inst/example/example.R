bb <- MutationTimeR:::refLengths[1:23]
bb$major_cn <- 3
bb$minor_cn <- 1
bb$clonal_frequency <- 0.8
bb$timing_param <- MutationTimeR:::defineMcnStates(bb, purity=0.8, clusters=data.frame(cluster=1, proportion=0.8, n_ssms=100))
for(i in seq_along(bb)){
	bb$timing_param[[i]][,"P.m.sX"] <- c(0.6,0.1,0.3)
	bb$timing_param[[i]][,"power.m.s"] <- rep(1,3)
}
bb$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6

vcf <- MutationTimeR:::simulateMutations(bb)
b <- bb[,c("major_cn","minor_cn","clonal_frequency")]

mt <- mutationTime(vcf, b, n.boot=10)
mcols(b) <- cbind(mcols(b),mt$T)

info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)

plotSample(vcf,b)
