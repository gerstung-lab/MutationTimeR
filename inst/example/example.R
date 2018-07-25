bb <- MutationTime.R:::refLengths[1:23]
bb$major_cn <- 3
bb$minor_cn <- 1
bb$clonal_frequency <- 0.8
bb$timing_param <- MutationTime.R:::defineMcnStates(bb, purity=0.8, clusters=data.frame(cluster=1, proportion=0.8, n_ssms=100))
for(i in seq_along(bb)){
	bb$timing_param[[i]][,"P.m.sX"] <- c(0.6,0.1,0.3)
	bb$timing_param[[i]][,"power.m.s"] <- rep(1,3)
}
bb$n.snv_mnv <- width(MutationTime.R:::refLengths[1:23])/1e6

v <- MutationTime.R:::simulateMutations(bb)
b <- bb[,c("major_cn","minor_cn","clonal_frequency")]

mt <- mutationTime(v, b, n.boot=10)
mcols(b) <- cbind(mcols(b),mt$T)

info(header(v)) <- rbind(info(header(v)),mtHeader())
info(v) <- cbind(info(v), mt$V)

plotSample(v,b)
