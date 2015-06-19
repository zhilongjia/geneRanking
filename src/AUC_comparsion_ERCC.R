#! /usr/bin/Rscript

library(pROC)
library(dplyr)
load("../results/runRanking.RData")
ERCC <- read.table("../data/ERCC_Controls_Analysis.txt", sep="\t", header=T, stringsAsFactors=FALSE)
ERCC <- ERCC[,c("ERCC.ID", "log2.Mix.1.Mix.2.")]
colnames(ERCC) <- c("ERCC.ID", "ERCC.logFC")
ERCC$ERCC.logFC <- -ERCC$ERCC.logFC

################################################################################
#AGR
gfold_AGR <- read.table("../results/AGR/gfold_AGR.rnk", sep="\t", row.names=1)
gfoldFC_AGR <- read.table("../results/AGR/gfoldFC_AGR.rnk", sep="\t", row.names=1)
pred_label_AGR <- data.frame(PoissonSeq=PoissonSeq_AGR$rnk[ERCC$ERCC.ID,], 
                             DNMF=DNMF_AGR$rnk[ERCC$ERCC.ID,],
                             edgeR=edgeR_AGR$rnk[ERCC$ERCC.ID,],
                             gfold=gfold_AGR[ERCC$ERCC.ID,], 
                             gfoldFC=gfoldFC_AGR[ERCC$ERCC.ID,], 
                             DESeq=DESeq_AGR$rnk[ERCC$ERCC.ID,])
#make NA as 0 since PoissonSeq filter small read counts internally. Other datasets 
#will be dealt with the same as AGR.
pred_label_AGR[which(is.na(pred_label_AGR$PoissonSeq)), 1]=0

################################################################################
#BGI
gfold_BGI <- read.table("../results/BGI/gfold_BGI.rnk", sep="\t", row.names=1)
gfoldFC_BGI <- read.table("../results/BGI/gfoldFC_BGI.rnk", sep="\t", row.names=1)
pred_label_BGI <- data.frame(PoissonSeq=PoissonSeq_BGI$rnk[ERCC$ERCC.ID,], 
                             DNMF=DNMF_BGI$rnk[ERCC$ERCC.ID,],
                             edgeR=edgeR_BGI$rnk[ERCC$ERCC.ID,],
                             gfold=gfold_BGI[ERCC$ERCC.ID,], 
                             gfoldFC=gfoldFC_BGI[ERCC$ERCC.ID,], 
                             DESeq=DESeq_BGI$rnk[ERCC$ERCC.ID,])
pred_label_BGI[which(is.na(pred_label_BGI$PoissonSeq)), 1]=0

################################################################################
#NWU
gfold_NWU <- read.table("../results/NWU/gfold_NWU.rnk", sep="\t", row.names=1)
gfoldFC_NWU <- read.table("../results/NWU/gfoldFC_NWU.rnk", sep="\t", row.names=1)
pred_label_NWU <- data.frame(PoissonSeq=PoissonSeq_NWU$rnk[ERCC$ERCC.ID,], 
                             DNMF=DNMF_NWU$rnk[ERCC$ERCC.ID,],
                             edgeR=edgeR_NWU$rnk[ERCC$ERCC.ID,],
                             gfold=gfold_NWU[ERCC$ERCC.ID,], 
                             gfoldFC=gfoldFC_NWU[ERCC$ERCC.ID,], 
                             DESeq=DESeq_NWU$rnk[ERCC$ERCC.ID,])
pred_label_NWU[which(is.na(pred_label_NWU$PoissonSeq)), 1]=0

################################################################################
#PSU
gfold_PSU <- read.table("../results/PSU/gfold_PSU.rnk", sep="\t", row.names=1)
gfoldFC_PSU <- read.table("../results/PSU/gfoldFC_PSU.rnk", sep="\t", row.names=1)
pred_label_PSU <- data.frame(PoissonSeq=PoissonSeq_PSU$rnk[ERCC$ERCC.ID,], 
                             DNMF=DNMF_PSU$rnk[ERCC$ERCC.ID,],
                             edgeR=edgeR_PSU$rnk[ERCC$ERCC.ID,],
                             gfold=gfold_PSU[ERCC$ERCC.ID,], 
                             gfoldFC=gfoldFC_PSU[ERCC$ERCC.ID,], 
                             DESeq=DESeq_PSU$rnk[ERCC$ERCC.ID,])
pred_label_PSU[which(is.na(pred_label_PSU$PoissonSeq)), 1]=0

################################################################################
#AUC with different cutoff
AUC_plot <- function(pred_label, ERCC, main) {
    cols <- c("red", "chocolate", "chartreuse", "deepskyblue", "darkgreen", "blue")
    auc_mat <- sapply(colnames(pred_label), 
                      function (x) {
                          sapply(c(0.58,1,2), function (y) auc(transmute(ERCC, label=ifelse(abs(ERCC.logFC)>=y, 1, 0))$label, abs(pred_label[,x])))
                      })
    plot(c(0.58,1,2), auc_mat[,1], type='n', main=main,
         xlab="logFC cutoff values", ylab="AUC", axes = FALSE,
         ylim=c(0.3,1))
    axis(side = 1, at = c(0.58, 1, 2))
    axis(side = 2, at = seq(0.2,1,0.1))
    for(i in seq(dim(auc_mat)[2])){
        lines(c(0.58,1,2), auc_mat[,i],
              lwd=3, col=cols[i], type="o")
    }
    legend("bottomleft", legend=colnames(auc_mat), col=cols, lwd=3,  cex=.6, ncol = 2, bty="n")
}
tiff(file="../results/Figure/AUC_ERCC.tiff", res=300, width = 5, height = 6, units="in", compression="lzw")
par(mfrow=c(2,2))
AUC_plot(pred_label_AGR, ERCC, main="AUCs for AGR (ERCC)")
AUC_plot(pred_label_BGI, ERCC, main="AUCs for BGI (ERCC)")
AUC_plot(pred_label_NWU, ERCC, main="AUCs for NWU (ERCC)")
AUC_plot(pred_label_PSU, ERCC, main="AUCs for PSU (ERCC)")
dev.off()
################################################################################
#Co-expression analysis
#AGR
# dat <- read.table("../data/AGR/AGR_readcount.txt", row.names=1)
# colnames(dat) = c(paste0("A", c(1,2,3,4)), paste0("B", c(1,2,3,4)), "GeneLength")
# dat <- as.matrix(dat[,-ncol(dat)])
# Sizefactors <- DESeq::estimateSizeFactorsForMatrix(dat)
# dat = sweep(dat, 2, Sizefactors, `/`)
# dat = dat[ERCC$ERCC.ID,]
# AGR_ERCC_rank <-  cbind(dat, apply(pred_label_AGR, 2, function(x) rank(x)))




