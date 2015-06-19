#! /usr/bin/Rscript

library(pROC)
library(dplyr)
load("../results/runRanking.RData")
taqman <- read.table("../data/seqc_taqman.txt", sep="\t", header=T)
taqman$Symbol <- as.character(taqman$Symbol)


################################################################################
#AGR
gfold_AGR <- read.table("../results/AGR/gfold_AGR.rnk", sep="\t", row.names=1)
gfoldFC_AGR <- read.table("../results/AGR/gfoldFC_AGR.rnk", sep="\t", row.names=1)
pred_label_AGR <- data.frame(PoissonSeq=PoissonSeq_AGR$rnk[taqman$Symbol,], 
                             DNMF=DNMF_AGR$rnk[taqman$Symbol,],
                             edgeR=edgeR_AGR$rnk[taqman$Symbol,],
                             gfold=gfold_AGR[taqman$Symbol,], 
                             gfoldFC=gfoldFC_AGR[taqman$Symbol,], 
                             DESeq=DESeq_AGR$rnk[taqman$Symbol,])

################################################################################
#BGI
gfold_BGI <- read.table("../results/BGI/gfold_BGI.rnk", sep="\t", row.names=1)
gfoldFC_BGI <- read.table("../results/BGI/gfoldFC_BGI.rnk", sep="\t", row.names=1)
pred_label_BGI <- data.frame(PoissonSeq=PoissonSeq_BGI$rnk[taqman$Symbol,], 
                             DNMF=DNMF_BGI$rnk[taqman$Symbol,],
                             edgeR=edgeR_BGI$rnk[taqman$Symbol,],
                             gfold=gfold_BGI[taqman$Symbol,], 
                             gfoldFC=gfoldFC_BGI[taqman$Symbol,], 
                             DESeq=DESeq_BGI$rnk[taqman$Symbol,])

################################################################################
#NWU
gfold_NWU <- read.table("../results/NWU/gfold_NWU.rnk", sep="\t", row.names=1)
gfoldFC_NWU <- read.table("../results/NWU/gfoldFC_NWU.rnk", sep="\t", row.names=1)
pred_label_NWU <- data.frame(PoissonSeq=PoissonSeq_NWU$rnk[taqman$Symbol,], 
                             DNMF=DNMF_NWU$rnk[taqman$Symbol,],
                             edgeR=edgeR_NWU$rnk[taqman$Symbol,],
                             gfold=gfold_NWU[taqman$Symbol,], 
                             gfoldFC=gfoldFC_NWU[taqman$Symbol,], 
                             DESeq=DESeq_NWU$rnk[taqman$Symbol,])

################################################################################
#PSU
gfold_PSU <- read.table("../results/PSU/gfold_PSU.rnk", sep="\t", row.names=1)
gfoldFC_PSU <- read.table("../results/PSU/gfoldFC_PSU.rnk", sep="\t", row.names=1)
pred_label_PSU <- data.frame(PoissonSeq=PoissonSeq_PSU$rnk[taqman$Symbol,], 
                             DNMF=DNMF_PSU$rnk[taqman$Symbol,],
                             edgeR=edgeR_PSU$rnk[taqman$Symbol,],
                             gfold=gfold_PSU[taqman$Symbol,], 
                             gfoldFC=gfoldFC_PSU[taqman$Symbol,], 
                             DESeq=DESeq_PSU$rnk[taqman$Symbol,])



################################################################################
#AUC with different cutoff
AUC_plot <- function(pred_label, taqman, main) {
    cols <- c("red", "chocolate", "chartreuse", "deepskyblue", "darkgreen", "blue")
    auc_mat <- sapply(colnames(pred_label), 
                      function (x) {
                          sapply(seq(0.5,2,0.1), function (y) auc(transmute(taqman, Symbol, label=ifelse(abs(logFC)>=y, 1, 0))$label, abs(pred_label[,x])))
                      })
    plot(seq(0.5,2,0.1), auc_mat[,1], type='n', main=main,
         xlab="logFC cutoff values", ylab="AUC",
         ylim=c(0.80,1))
    
    for(i in seq(dim(auc_mat)[2])){
        lines(seq(0.5,2,0.1), auc_mat[,i],
              lwd=3, col=cols[i])
    }
    legend("bottomright", legend=colnames(auc_mat), col=cols, lwd=3,  cex=.6, ncol = 2, bty="n")
}
tiff(file="../results/Figure/AUC_taqman.tiff", res=300, width = 5, height = 6, units="in", compression="lzw")
par(mfrow=c(2,2))
AUC_plot(pred_label_AGR, taqman, main="AUCs for AGR")
AUC_plot(pred_label_BGI, taqman, main="AUCs for BGI")
AUC_plot(pred_label_NWU, taqman, main="AUCs for NWU")
AUC_plot(pred_label_PSU, taqman, main="AUCs for PSU")
dev.off()
