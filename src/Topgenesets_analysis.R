#! /usr/bin/Rscript

library(reshape2)
library(dplyr)

geneset <- read.table("../results/topgenesets.xls", sep="\t", header=TRUE, stringsAsFactors=FALSE)


#c2
venn <- list()
for (i in c("DESeq", "DNMF", "edgeR", "gfold", "gfoldFC", "PossionSeq")){
    j <- paste( i, "(C2CP)")
    print (j)
    venn[[j]] <- list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c2cp", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c2cp", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c2cp", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c2cp", method==i)[,4:13]) )
}



#c2kegg
# venn <- list()
for (i in c("DESeq", "DNMF", "edgeR", "gfold", "gfoldFC", "PossionSeq")){
    j <- paste( i, "(C2KEGG)")
    print (j)
    venn[[j]] <- venn.diagram(list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c2kegg", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c2kegg", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c2kegg", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c2kegg", method==i)[,4:13]) ), 
                              main=j, height = 300, width = 300, filename=NULL,  sub.cex=0.5);
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
    j <- paste( i, "(C5BP)")
    print (j)
    venn[[j]] <- venn.diagram(list(AGR=t(dplyr::filter(geneset, dataset=="AGR", geneset=="c5bp", method==i)[,4:13]),
                                   BGI=t(dplyr::filter(geneset, dataset=="BGI", geneset=="c5bp", method==i)[,4:13]),
                                   NWU=t(dplyr::filter(geneset, dataset=="NWU", geneset=="c5bp", method==i)[,4:13]),
                                   PSU=t(dplyr::filter(geneset, dataset=="PSU", geneset=="c5bp", method==i)[,4:13]) ), 
                              main=j, height = 300, width = 300, filename=NULL,  sub.cex=0.5);
}

# grid.arrange(gTree(children=venn[[1]]),
#              gTree(children=venn[[2]]),
#              gTree(children=venn[[3]]),
#              gTree(children=venn[[4]]),
#              gTree(children=venn[[5]]),
#              gTree(children=venn[[6]]),
#              ncol=1, nrow=6)

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

##################################################################################################
#The overlapped gene sets by DNMF. top10

C2CP <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c2cp", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C2CP), file="../results/DNMF_C2CP_intersect.txt", row.names=F, col.names=F, quote=F)

C2KEGG <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c2kegg", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C2KEGG), file="../results/DNMF_C2KEGG_intersect.txt", row.names=F, col.names=F, quote=F)

C5 <- as.list(as.data.frame(t(dplyr::filter(geneset, geneset=="c5bp", method=="DNMF")[4:13])))
write.table(Reduce(intersect, C5), file="../results/DNMF_C5_intersect.txt", row.names=F, col.names=F, quote=F)




################################################################################
#c2cp
c2cp <- dplyr::filter(geneset, geneset=="c2cp") %>% melt(1:3) %>% count(value) %>%  dplyr::filter(n==24)

AGRc2cp <- dplyr::filter(geneset, dataset=="AGR", geneset=="c2cp") %>% reshape2::melt(1:3)
AGRc2cp <- dplyr::filter(AGRc2cp, value %in% c2cp$value)[,3:5] %>%  reshape2::dcast(value~method, value.var="variable")

BGIc2cp <- dplyr::filter(geneset, dataset=="BGI", geneset=="c2cp") %>% reshape2::melt(1:3)
BGIc2cp <- dplyr::filter(BGIc2cp, value %in% c2cp$value)[,3:5] %>%  reshape2::dcast(value~method, value.var="variable")


NWUc2cp <- dplyr::filter(geneset, dataset=="NWU", geneset=="c2cp") %>% reshape2::melt(1:3)
NWUc2cp <- dplyr::filter(NWUc2cp, value %in% c2cp$value)[,3:5] %>%  reshape2::dcast(value~method, value.var="variable")

PSUc2cp <- dplyr::filter(geneset, dataset=="PSU", geneset=="c2cp") %>% reshape2::melt(1:3)
PSUc2cp <- dplyr::filter(PSUc2cp, value %in% c2cp$value)[,3:5] %>%  reshape2::dcast(value~method, value.var="variable")


