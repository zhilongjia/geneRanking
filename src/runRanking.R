#! /usr/bin/Rscript

source("./geneRanking.R")
################################################################################
#AGR
dat <- read.table("../data/AGR/AGR_readcount.txt", row.names=1)
colnames(dat) = c(paste0("A", c(1,2,3,4)), paste0("B", c(1,2,3,4)), "GeneLength")
dat <- as.matrix(dat[,-ncol(dat)])
samplelabel = as.numeric(rep(seq(1:2), each=ncol(dat)/2))

DESeq_AGR <- runDESeq(dat, samplelabel)
DNMF_AGR <- runDNMF(dat, samplelabel)
edgeR_AGR <- runedgeR(dat, samplelabel)
PoissonSeq_AGR <- runPoissonSeq(dat, samplelabel)

write.table(DESeq_AGR$rnk, "../results/AGR/DESeq_AGR.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(DNMF_AGR$rnk, "../results/AGR/DNMF_AGR.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(edgeR_AGR$rnk, "../results/AGR/edgeR_AGR.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(PoissonSeq_AGR$rnk, "../results/AGR/PoissonSeq_AGR.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
rm(dat, samplelabel)

################################################################################
#NWU
dat <- read.table("../data/NWU/NWU_readcount.txt", row.names=1)
colnames(dat) = c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
dat <- as.matrix(dat[,-ncol(dat)])
samplelabel = as.numeric(rep(seq(1:2), each=ncol(dat)/2))

strt<-Sys.time()
DESeq_NWU <- runDESeq(dat, samplelabel)
print(Sys.time()-strt)
strt<-Sys.time()
DNMF_NWU <- runDNMF(dat, samplelabel)
print(Sys.time()-strt)
strt<-Sys.time()
edgeR_NWU <- runedgeR(dat, samplelabel)
print(Sys.time()-strt)
strt<-Sys.time()
PoissonSeq_NWU <- runPoissonSeq(dat, samplelabel)
print(Sys.time()-strt)


write.table(DESeq_NWU$rnk, "../results/NWU/DESeq_NWU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(DNMF_NWU$rnk, "../results/NWU/DNMF_NWU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(edgeR_NWU$rnk, "../results/NWU/edgeR_NWU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(PoissonSeq_NWU$rnk, "../results/NWU/PoissonSeq_NWU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
rm(dat, samplelabel)

#NWU_raw (only testing time consumption)
dat <- read.table("../data/NWU_raw/NWU_raw_readcount.txt", sep="\t", row.names=1)
dat <- as.matrix(dat[,-ncol(dat)])
samplelabel = as.numeric(rep(seq(1:2), each=ncol(dat)/2))
strt<-Sys.time()
DNMF_NWU_raw <- runDNMF(dat, samplelabel)
print(Sys.time()-strt)
rm(dat, samplelabel)
rnk(dnmf_result, "../results/NWU_raw/DNMF_NWU_raw.rnk")

################################################################################
#BGI
dat <- read.table("../data/BGI/BGI_readcount.txt", row.names=1)
colnames(dat) = c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
dat <- as.matrix(dat[,-ncol(dat)])
samplelabel = as.numeric(rep(seq(1:2), each=ncol(dat)/2))

DESeq_BGI <- runDESeq(dat, samplelabel)
DNMF_BGI <- runDNMF(dat, samplelabel)
edgeR_BGI <- runedgeR(dat, samplelabel)
PoissonSeq_BGI <- runPoissonSeq(dat, samplelabel)
write.table(DESeq_BGI$rnk, "../results/BGI/DESeq_BGI.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(DNMF_BGI$rnk, "../results/BGI/DNMF_BGI.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(edgeR_BGI$rnk, "../results/BGI/edgeR_BGI.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(PoissonSeq_BGI$rnk, "../results/BGI/PoissonSeq_BGI.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
rm(dat, samplelabel)

################################################################################
#PSU
dat <- read.table("../data/PSU/PSU_readcount.txt", row.names=1)
colnames(dat) = c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
dat <- as.matrix(dat[,-ncol(dat)])
samplelabel = as.numeric(rep(seq(1:2), each=ncol(dat)/2))

DESeq_PSU <- runDESeq(dat, samplelabel)
DNMF_PSU <- runDNMF(dat, samplelabel)
edgeR_PSU <- runedgeR(dat, samplelabel)
PoissonSeq_PSU <- runPoissonSeq(dat, samplelabel)
write.table(DESeq_PSU$rnk, "../results/PSU/DESeq_PSU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(DNMF_PSU$rnk, "../results/PSU/DNMF_PSU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(edgeR_PSU$rnk, "../results/PSU/edgeR_PSU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
write.table(PoissonSeq_PSU$rnk, "../results/PSU/PoissonSeq_PSU.rnk", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")
rm(dat, samplelabel)

save.image("../results/runRanking.RData")

