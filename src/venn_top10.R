#! /usr/bin/Rscript

#Plotting the venn Diagram of four datasets for c2cp, c2kegg and c5 respectively.

library(VennDiagram)
library(gridExtra)

geneset <- read.table("../results/topgenesets.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)

geneset$dataset <- rep(c("AGR", "BGI", "NWU", "PSU"), each=18)
geneset$geneset <- rep(rep(c("c2cp", "c2kegg", "c5bp"), each=6), 4)
geneset$method <- rep(c("DESeq", "DNMF", "edgeR", "gfoldFC", "gfold", "PossionSeq"), 12)
geneset <- geneset[,c("dataset", "geneset", "method", paste0("V", 2:11))]
write.table(geneset, file="../results/topgenesets.xls", sep="\t", quote=FALSE, row.names = FALSE)


#c2
venn <- list()
for (i in c("DESeq", "DNMF", "edgeR", "gfold", "gfoldFC", "PossionSeq")){
    j <- paste( i, "(C2CP)")
    print (j)
    venn[[j]] <- venn.diagram(list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c2cp", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c2cp", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c2cp", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c2cp", method==i)[,4:13]) ), 
                              main=j, height = 300, width = 300, filename=NULL,  main.cex=0.9, cex=0.7, cat.cex=0.7);
}

# grid.arrange(gTree(children=venn[[1]]),
#              gTree(children=venn[[2]]),
#              gTree(children=venn[[3]]),
#              gTree(children=venn[[4]]),
#              gTree(children=venn[[5]]),
#              gTree(children=venn[[6]]),
#              ncol=1, nrow=6)

#c2kegg
# venn <- list()
for (i in c("DESeq", "DNMF", "edgeR", "gfold", "gfoldFC", "PossionSeq")){
    j <- paste( i, "(C2KEGG)")
    print (j)
    venn[[j]] <- venn.diagram(list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c2kegg", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c2kegg", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c2kegg", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c2kegg", method==i)[,4:13]) ), 
                              main=j, height = 300, width = 300, filename=NULL,  main.cex=0.9, cex=0.7, cat.cex=0.7);
}

# grid.arrange(gTree(children=venn[[1]]),
#              gTree(children=venn[[2]]),
#              gTree(children=venn[[3]]),
#              gTree(children=venn[[4]]),
#              gTree(children=venn[[5]]),
#              gTree(children=venn[[6]]),
#              ncol=1, nrow=6)

#c5bp
# venn <- list()
for (i in c("DESeq", "DNMF", "edgeR", "gfold", "gfoldFC", "PossionSeq")){
    j <- paste( i, "(C5)")
    print (j)
    venn[[j]] <- venn.diagram(list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c5bp", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c5bp", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c5bp", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c5bp", method==i)[,4:13]) ), 
                              main=j, height = 300, width = 300, filename=NULL,  main.cex=0.9, cex=0.7, cat.cex=0.7);
}

# grid.arrange(gTree(children=venn[[1]]),
#              gTree(children=venn[[2]]),
#              gTree(children=venn[[3]]),
#              gTree(children=venn[[4]]),
#              gTree(children=venn[[5]]),
#              gTree(children=venn[[6]]),
#              ncol=1, nrow=6)

tiff(file="../results/Figure/Venn_top10.tiff", res=300, width = 6, height = 7, units="in", compression="lzw")
grid.arrange(gTree(children=venn[[1]]),
             gTree(children=venn[[7]]),
             gTree(children=venn[[13]]),
             gTree(children=venn[[2]]),
             gTree(children=venn[[8]]),
             gTree(children=venn[[14]]),
             gTree(children=venn[[3]]),
             gTree(children=venn[[9]]),
             gTree(children=venn[[15]]),
             gTree(children=venn[[4]]),
             gTree(children=venn[[10]]),
             gTree(children=venn[[16]]),
             gTree(children=venn[[5]]),
             gTree(children=venn[[11]]),
             gTree(children=venn[[17]]),
             gTree(children=venn[[6]]),
             gTree(children=venn[[12]]),
             gTree(children=venn[[18]]),
             ncol=3, nrow=6)
dev.off()

##################################################################################################
#The overlapped gene sets by DNMF. top10

C2CP <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c2cp", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C2CP), file="../results/DNMF_C2CP_intersect.txt", row.names=F, col.names=F, quote=F)

C2KEGG <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c2kegg", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C2KEGG), file="../results/DNMF_C2KEGG_intersect.txt", row.names=F, col.names=F, quote=F)

C5 <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c5bp", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C5), file="../results/DNMF_C5_intersect.txt", row.names=F, col.names=F, quote=F)


