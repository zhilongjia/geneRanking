#! /usr/bin/Rscript

#DESeq
DESeq_LE <- read.table("../results/LeadingEdge_AGR/DESeq_C2CP_leading_edge_matrix_for_results.8.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
DESeq_LEG <- names(which(colSums(DESeq_LE)==max(colSums(DESeq_LE))))

#DNMF
DNMF_LE <- read.table("../results/LeadingEdge_AGR/DNMF_C2CP_leading_edge_matrix_for_results.11.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
DNMF_LEG <- names(which(colSums(DNMF_LE)==max(colSums(DNMF_LE))))

#edgeR
edgeR_LE <- read.table("../results/LeadingEdge_AGR/edgeR_C2CP_leading_edge_matrix_for_results.13.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
edgeR_LEG <- names(which(colSums(edgeR_LE)==max(colSums(edgeR_LE))))

#gfold
gfold_LE <- read.table("../results/LeadingEdge_AGR/gfold_C2CP_leading_edge_matrix_for_results.15.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
gfold_LEG <- names(which(colSums(gfold_LE)==max(colSums(gfold_LE))))

#gfoldFC
gfoldFC_LE <- read.table("../results/LeadingEdge_AGR/gfoldFC_C2CP_leading_edge_matrix_for_results.14.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
gfoldFC_LEG <- names(which(colSums(gfoldFC_LE)==max(colSums(gfoldFC_LE))))

#PossionSeq
PossionSeq_LE <- read.table("../results/LeadingEdge_AGR/PossionSeq_C2CP_leading_edge_matrix_for_results.12.gct", sep="\t", skip=2, header=TRUE)[,-c(1,2)]
PossionSeq_LEG <- names(which(colSums(PossionSeq_LE)==max(colSums(PossionSeq_LE))))

LEG <- list(DESeq_LEG=DESeq_LEG, DNMF_LEG=DNMF_LEG, edgeR_LEG=edgeR_LEG, gfold_LEG=gfold_LEG, gfoldFC_LEG=gfoldFC_LEG, PossionSeq_LEG=PossionSeq_LEG)
write.table(unlist(LEG), file="../results/LEG.txt", sep="\t", col.names = F, quote=F)

#Venn Diagram
vLEG <- list(DESeq=DESeq_LEG, "DNMF/gfold"=DNMF_LEG, edgeR=edgeR_LEG, "gfoldFC/PossionSeq"=gfoldFC_LEG)

#Venn Diagram
VennDiagram::venn.diagram(vLEG, filename="2.tiff", res=300, width = 6, height = 6, units="in", compression="lzw")


save.image("../results/LeadingEdgeAnalysis.RData")
