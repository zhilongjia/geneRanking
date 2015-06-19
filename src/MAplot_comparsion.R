#! /usr/bin/Rscript

library(limma)
library(dplyr)
load("../results/runRanking.RData")

################################################################################
#AGR
data <- read.table("../data/AGR/AGR_readcount.txt", sep="\t", row.names=1)
data <- as.matrix(data[,-ncol(data)])
gfold_AGR <- read.table("../results/AGR/gfold_AGR.rnk", sep="\t", row.names=1)
gfoldFC_AGR <- read.table("../results/AGR/gfoldFC_AGR.rnk", sep="\t", row.names=1)
colnames(gfoldFC_AGR) = colnames(gfold_AGR) <- "rnk"
rnkList_AGR <- list(PoissonSeq=PoissonSeq_AGR$rnk, 
                    DNMF=DNMF_AGR$rnk, 
                    edgeR=edgeR_AGR$rnk, 
                    gfold=gfold_AGR, 
                    gfoldFC=gfoldFC_AGR, 
                    DESeq=DESeq_AGR$rnk)

#Plot MA (top 1000 genes)
tiff(file="../results/Figure/MA_AGR.tiff", res=300, width = 5.2, height = 4, units="in", compression="lzw")
par(mfrow=c(2,3))
for (i in 1: length(rnkList_AGR)) {
    rnk=rnkList_AGR[[i]]
    rnk$Symbol <- rownames(rnk)
    gene <- as.character(top_n(rnk, 1000, abs(rnk))$Symbol)
    dataMean <- as.data.frame(data[as.character(rnk$Symbol),]+1)
    dataMean$Symbol <- rownames(dataMean)
    # dataMean <- dplyr::transmute(dataMean, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    dataMean <- dplyr::transmute(dataMean, Symbol, meanA=(V2+V3+V4+V5)/4, meanB=(V6+V7+V8+V9)/4)
    
    status <- ifelse(as.character(rnk$Symbol) %in% gene, "1", "0")
    MA <- new("MAList")
    MA$A <- 1/2 * (log2(dataMean[,"meanA"]) + log2(dataMean[,"meanB"]))
    MA$M <- log2(dataMean[,"meanB"]) - log2(dataMean[,"meanA"])
    limma::plotMA(MA, main=names(rnkList_AGR)[i], status=status, legend=FALSE)
}
dev.off()
rm(data)

################################################################################
#NWU
data <- read.table("../data/NWU/NWU_readcount.txt", sep="\t", row.names=1)
data <- as.matrix(data[,-ncol(data)])
gfold_NWU <- read.table("../results/NWU/gfold_NWU.rnk", sep="\t", row.names=1)
gfoldFC_NWU <- read.table("../results/NWU/gfoldFC_NWU.rnk", sep="\t", row.names=1)
colnames(gfoldFC_NWU) = colnames(gfold_NWU) <- "rnk"
rnkList_NWU <- list(PoissonSeq=PoissonSeq_NWU$rnk, 
                    DNMF=DNMF_NWU$rnk, 
                    edgeR=edgeR_NWU$rnk, 
                    gfold=gfold_NWU, 
                    gfoldFC=gfoldFC_NWU, 
                    DESeq=DESeq_NWU$rnk)

#Plot MA (top 1000 genes)
tiff(file="../results/Figure/MA_NWU.tiff", res=300, width = 5.2, height = 4, units="in", compression="lzw")
par(mfrow=c(2,3))
for (i in 1: length(rnkList_NWU)) {
    rnk=rnkList_NWU[[i]]
    rnk$Symbol <- rownames(rnk)
    gene <- as.character(top_n(rnk, 1000, abs(rnk))$Symbol)
    dataMean <- as.data.frame(data[as.character(rnk$Symbol),]+1)
    dataMean$Symbol <- rownames(dataMean)
    # dataMean <- dplyr::transmute(dataMean, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    dataMean <- dplyr::transmute(dataMean, Symbol, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    
    status <- ifelse(as.character(rnk$Symbol) %in% gene, "1", "0")
    MA <- new("MAList")
    MA$A <- 1/2 * (log2(dataMean[,"meanA"]) + log2(dataMean[,"meanB"]))
    MA$M <- log2(dataMean[,"meanB"]) - log2(dataMean[,"meanA"])
    limma::plotMA(MA, main=names(rnkList_NWU)[i], status=status, legend=FALSE)
}
dev.off()
rm(data)
################################################################################
#BGI
data <- read.table("../data/BGI/BGI_readcount.txt", sep="\t", row.names=1)
data <- as.matrix(data[,-ncol(data)])
gfold_BGI <- read.table("../results/BGI/gfold_BGI.rnk", sep="\t", row.names=1)
gfoldFC_BGI <- read.table("../results/BGI/gfoldFC_BGI.rnk", sep="\t", row.names=1)
colnames(gfoldFC_BGI) = colnames(gfold_BGI) <- "rnk"
rnkList_BGI <- list(PoissonSeq=PoissonSeq_BGI$rnk, 
                    DNMF=DNMF_BGI$rnk, 
                    edgeR=edgeR_BGI$rnk, 
                    gfold=gfold_BGI, 
                    gfoldFC=gfoldFC_BGI, 
                    DESeq=DESeq_BGI$rnk)

#Plot MA (top 1000 genes)
tiff(file="../results/Figure/MA_BGI.tiff", res=300, width = 5.2, height = 4, units="in", compression="lzw")
par(mfrow=c(2,3))
for (i in 1: length(rnkList_BGI)) {
    rnk=rnkList_BGI[[i]]
    rnk$Symbol <- rownames(rnk)
    gene <- as.character(top_n(rnk, 1000, abs(rnk))$Symbol)
    dataMean <- as.data.frame(data[as.character(rnk$Symbol),]+1)
    dataMean$Symbol <- rownames(dataMean)
    # dataMean <- dplyr::transmute(dataMean, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    dataMean <- dplyr::transmute(dataMean, Symbol, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    
    status <- ifelse(as.character(rnk$Symbol) %in% gene, "1", "0")
    MA <- new("MAList")
    MA$A <- 1/2 * (log2(dataMean[,"meanA"]) + log2(dataMean[,"meanB"]))
    MA$M <- log2(dataMean[,"meanB"]) - log2(dataMean[,"meanA"])
    limma::plotMA(MA, main=names(rnkList_BGI)[i], status=status, legend=FALSE)
}
dev.off()
rm(data)
################################################################################
#PSU
data <- read.table("../data/PSU/PSU_readcount.txt", sep="\t", row.names=1)
data <- as.matrix(data[,-ncol(data)])
gfold_PSU <- read.table("../results/PSU/gfold_PSU.rnk", sep="\t", row.names=1)
gfoldFC_PSU <- read.table("../results/PSU/gfoldFC_PSU.rnk", sep="\t", row.names=1)
colnames(gfoldFC_PSU) = colnames(gfold_PSU) <- "rnk"
rnkList_PSU <- list(PoissonSeq=PoissonSeq_PSU$rnk, 
                    DNMF=DNMF_PSU$rnk, 
                    edgeR=edgeR_PSU$rnk, 
                    gfold=gfold_PSU, 
                    gfoldFC=gfoldFC_PSU, 
                    DESeq=DESeq_PSU$rnk)

#Plot MA (top 1000 genes)
tiff(file="../results/Figure/MA_PSU.tiff", res=300, width = 5.2, height = 4, units="in", compression="lzw")
par(mfrow=c(2,3))
for (i in 1: length(rnkList_PSU)) {
    rnk=rnkList_PSU[[i]]
    rnk$Symbol <- rownames(rnk)
    gene <- as.character(top_n(rnk, 1000, abs(rnk))$Symbol)
    dataMean <- as.data.frame(data[as.character(rnk$Symbol),]+1)
    dataMean$Symbol <- rownames(dataMean)
    # dataMean <- dplyr::transmute(dataMean, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    dataMean <- dplyr::transmute(dataMean, Symbol, meanA=(V2+V3+V4+V5+V6)/5, meanB=(V7+V8+V9+V10+V11)/5)
    
    status <- ifelse(as.character(rnk$Symbol) %in% gene, "1", "0")
    MA <- new("MAList")
    MA$A <- 1/2 * (log2(dataMean[,"meanA"]) + log2(dataMean[,"meanB"]))
    MA$M <- log2(dataMean[,"meanB"]) - log2(dataMean[,"meanA"])
    limma::plotMA(MA, main=names(rnkList_PSU)[i], status=status, legend=FALSE)
}
dev.off()
