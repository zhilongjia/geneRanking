#! /usr/bin/R

library(seqc)
ls(2)

data(ILM_refseq_gene_AGR)
sampleNames <- colnames(ILM_refseq_gene_AGR)
AGR <- cbind(rowSums(ILM_refseq_gene_AGR[,which(grepl("A_1", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("A_2", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("A_3", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("A_4", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("B_1", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("B_2", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("B_3", sampleNames))]),
      rowSums(ILM_refseq_gene_AGR[,which(grepl("B_4", sampleNames))]),
      ILM_refseq_gene_AGR$GeneLength)

colnames(AGR) <- c(paste0("A", c(1,2,3,4)), paste0("B", c(1,2,3,4)), "GeneLength")
rownames(AGR) <- ILM_refseq_gene_AGR$Symbol
AGR <- AGR[!is.na(rownames(AGR)),]
#delete two KIR3DL3
AGR <- AGR[which(rownames(AGR)!="KIR3DL3"),]

################################################################################
#BGI
data(ILM_refseq_gene_BGI)
sampleNames <- colnames(ILM_refseq_gene_BGI)
BGI <- cbind(rowSums(ILM_refseq_gene_BGI[,which(grepl("A_1", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("A_2", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("A_3", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("A_4", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("A_5", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("B_1", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("B_2", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("B_3", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("B_4", sampleNames))]),
             rowSums(ILM_refseq_gene_BGI[,which(grepl("B_5", sampleNames))]),
             ILM_refseq_gene_BGI$GeneLength)

colnames(BGI) <- c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
rownames(BGI) <- ILM_refseq_gene_BGI$Symbol
BGI <- BGI[!is.na(rownames(BGI)),]
#delete two KIR3DL3
which(rownames(BGI)=="KIR3DL3")
BGI <- BGI[which(rownames(BGI)!="KIR3DL3"),]

################################################################################
#NWU
data(LIF_refseq_gene_NWU)
sampleNames <- colnames(LIF_refseq_gene_NWU)
NWU <- cbind(rowSums(LIF_refseq_gene_NWU[,which(grepl("A_1", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("A_2", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("A_3", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("A_4", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("A_5", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("B_1", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("B_2", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("B_3", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("B_4", sampleNames))]),
             rowSums(LIF_refseq_gene_NWU[,which(grepl("B_5", sampleNames))]),
             LIF_refseq_gene_NWU$GeneLength)

colnames(NWU) <- c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
rownames(NWU) <- LIF_refseq_gene_NWU$Symbol
NWU <- NWU[!is.na(rownames(NWU)),]
#delete two KIR3DL3
NWU <- NWU[which(rownames(NWU)!="KIR3DL3"),]

################################################################################
#PSU
data(LIF_refseq_gene_PSU)
sampleNames <- colnames(LIF_refseq_gene_PSU)
PSU <- cbind(rowSums(LIF_refseq_gene_PSU[,which(grepl("A_1", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("A_2", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("A_3", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("A_4", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("A_5", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("B_1", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("B_2", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("B_3", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("B_4", sampleNames))]),
             rowSums(LIF_refseq_gene_PSU[,which(grepl("B_5", sampleNames))]),
             LIF_refseq_gene_PSU$GeneLength)

colnames(PSU) <- c(paste0("A", c(1,2,3,4,5)), paste0("B", c(1,2,3,4,5)), "GeneLength")
rownames(PSU) <- LIF_refseq_gene_PSU$Symbol
PSU <- PSU[!is.na(rownames(PSU)),]
#delete two KIR3DL3
PSU <- PSU[which(rownames(PSU)!="KIR3DL3"),]


write.table(AGR, file="../data/AGR/AGR_readcount.txt", sep="\t", col.names=F, quote=F)
write.table(BGI, file="../data/BGI/BGI_readcount.txt", sep="\t", col.names=F, quote=F)
write.table(NWU, file="../data/NWU/NWU_readcount.txt", sep="\t", col.names=F, quote=F)
write.table(PSU, file="../data/PSU/PSU_readcount.txt", sep="\t", col.names=F, quote=F)

################################################################################
#qPCR data
data(taqman)
taqmanData <- taqman[, grep("A|B", colnames(taqman))]
taqmanData$Symbol <- taqman$Symbol

taqmanData <- dplyr::transmute(taqmanData, Symbol, meanA=(A.A1_value + A.A2_value + A.A3_value + A.A4_value)/4, meanB=(B.B1_value + B.B2_value + B.B3_value + B.B4_value)/4) %>% 
    dplyr::mutate(logFC=log2(meanB/meanA))

write.table(taqmanData, file="../data/seqc_taqman.txt", sep="\t", row.names=F, quote=F)

################################################################################
#NWU_raw
NWU_raw <- as.matrix(cbind(LIF_refseq_gene_NWU[,grep("^A|B", sampleNames)], LIF_refseq_gene_NWU$GeneLength))

colnames(NWU_raw) <- c(paste0("A", 1:60), paste0("B", 1:60), "GeneLength")
rownames(NWU_raw) <- LIF_refseq_gene_NWU$Symbol
NWU_raw <- NWU_raw[!is.na(rownames(NWU_raw)),]
#delete two KIR3DL3
NWU_raw <- NWU_raw[which(rownames(NWU_raw)!="KIR3DL3"),]

write.table(NWU_raw, file="../data/NWU_raw/NWU_raw_readcount.txt", sep="\t", col.names=F, quote=F)

