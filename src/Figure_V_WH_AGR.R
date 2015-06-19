#! /usr/bin/Rscript
#plot the Fig 1

require(pheatmap)


###################################################
load("../results/runRanking.RData")
dnmf_result <- DNMF_AGR$res
data <- dnmf_result$V
# colnames(data) = colnames(dnmf_result$H) = c(paste0("A", c(1,2,3,4)), paste0("B", c(1,2,3,4)))

color = colorRampPalette(c("navy", "white", "firebrick3"))(50)

# data <- t(apply(data, 1, function(x){x/max(x)}))
data <- apply(data, 2, function(x) {x/max(x)})
#data <- t(scale(t(data), center = FALSE, scale = TRUE))
tiff(file="../results/Figure/V.tiff", res=600, width = 40, height = 80, units="mm")
pheatmap(data, color=color, scale="none", cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, cellwidth=10, main="V",legend=F)
dev.off()


tiff(file="../results/Figure/V_legend.tiff", res=600, width = 50, height = 80, units="mm")
pheatmap(data, color=color, scale="none", cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, cellwidth=10, main="V")
dev.off()
################################################################################
W <- dnmf_result$W
W <- apply(W, 2, function(x) {x/max(x)})
# W <- t(apply(W, 1, function(x){x/max(x)}))
#W <- t(scale(t(W), center = FALSE, scale = TRUE))
colnames(W) <- c("Down", "Up")

tiff(file="../results/Figure/W.tiff", res=600, width = 15, height = 80, units="mm")
pheatmap(W, color=color, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, cellwidth=10, main="W",legend=F)
dev.off()
################################################################################
H <- dnmf_result$H
H <- apply(H, 2, function(x) {x/max(x)})

# H <- t(apply(H, 1, function(x){x/max(x)}))
#H <- t(scale(t(H), center = FALSE, scale = TRUE))
rownames(H) <- c("Down", "Up")

tiff(file="../results/Figure/H.tiff", res=600, width = 40, height = 30, units="mm")
pheatmap(H, color=color, show_rownames=TRUE, cluster_rows=FALSE, cluster_cols=FALSE,  cellheight=15, cellwidth=8, main="H", legend=F)
dev.off()

################################################################################
#Classifing AGR samples with top 10 ranked genes 
rnk <- dnmf_result$rnk
gn <- names(sort(abs(rnk), decreasing = TRUE)[1:10])
tiff(file="../results/Figure/Heatmap_wedV.tiff", res=600, width = 60, height = 80, units="mm", compression="lzw")
pheatmap(data[gn,], color=color, scale="none", cluster_rows=FALSE, cluster_cols=TRUE, show_rownames=TRUE, cellwidth=10)
dev.off()




