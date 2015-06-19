runDESeq <- function(data, condition){
    
    library(DESeq)
    condition <- as.factor(condition)
    cds = newCountDataSet(data, condition)
    
    #normalisation and Variance estimation
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions( cds )
    #significant genes
    strt<-Sys.time()
    res = nbinomTest (cds, 1, 2)
    print(Sys.time()-strt)
    
    #using -log2(padj) with the sign of log2FoldChange
    res$pval[which(res$pval==0)] = min(res$pval[which(res$pval!=0)],na.rm=T)
    res$rnk <- - sign(res$log2FoldChange) * log2(res$pval)
    res[which(is.na(res$rnk)),"rnk"] <- median(res$rnk, na.rm=TRUE)
    logFC = ifelse (res[which(res$rnk == max(res$rnk)),"log2FoldChange"] == Inf, 
                    max(res[which(res$rnk == max(res$rnk)),"log2FoldChange"][which(is.finite(res[which(res$rnk == max(res$rnk)),"log2FoldChange"]))]),
                    res[which(res$rnk == max(res$rnk)),"log2FoldChange"])
    res[which(res$rnk == max(res$rnk)),"rnk"] <- res[which(res$rnk == max(res$rnk)),"rnk"] + logFC
    
    logFC = ifelse (res[which(res$rnk == min(res$rnk)),"log2FoldChange"] == -Inf, 
                    min(res[which(res$rnk == min(res$rnk)),"log2FoldChange"][which(is.finite(res[which(res$rnk == min(res$rnk)),"log2FoldChange"]))]),
                    res[which(res$rnk == min(res$rnk)),"log2FoldChange"])
    res[which(res$rnk == min(res$rnk)),"rnk"] <- res[which(res$rnk == min(res$rnk)),"rnk"] + logFC
    
    rnk <- res[, c("id", "rnk")]
    rnk[which(is.nan(rnk$rnk)),"rnk"] = 0
    rownames(rnk) <- rnk$id
    rnk=rnk[,c("rnk"), drop=FALSE]
    return (list(rnk=rnk, res=res))
}

################################################################################
runDNMF <- function(data, trainlabel, r = 2){
    library(DNMF)
    library(DESeq)
    
    ###################################
    #normalization
    Sizefactors <- estimateSizeFactorsForMatrix(data)
    data = sweep(data, 2, Sizefactors, `/`)
    
    ##########################
    # perform DNMF, like V=WH 
    dnmf_result = DNMF(data, trainlabel, r)
    rnk = as.data.frame(dnmf_result$rnk)
    colnames(rnk) <- "rnk"
    return (list(rnk=rnk, res=dnmf_result))
}

################################################################################
runedgeR <- function(data, samplelabel){
    library(edgeR)
    
    data <- DGEList(counts = data, group=samplelabel)
    data <- calcNormFactors(data)
    
    data <- estimateCommonDisp(data)
    data <- estimateTagwiseDisp(data)
    et <- exactTest(data)
    
    #generate the rnk file for GSEA
    etg <- et$table
    etg$PValue[which(etg$PValue==0)] = min(etg$PValue[which(etg$PValue!=0)],na.rm=T)
    etg$rnk <- - sign(etg$logFC) * log2(etg$PValue)
    etg[which(etg$rnk == max(etg$rnk)),"rnk"] <- etg[which(etg$rnk == max(etg$rnk)),"rnk"] + etg[which(etg$rnk == max(etg$rnk)),"logFC"]
    etg[which(etg$rnk == min(etg$rnk)),"rnk"] <- etg[which(etg$rnk == min(etg$rnk)),"rnk"] + etg[which(etg$rnk == min(etg$rnk)),"logFC"]
    
    rnk <- etg[, "rnk", drop=FALSE]
    return (list(rnk=rnk, res=etg))
}

################################################################################
runPoissonSeq <- function(count.dat, conditions){
    library("PoissonSeq")

    dat <- list(n=count.dat,
                y=conditions,
                type='twoclass',
                pair=FALSE,
                gname=rownames(count.dat))
    para <- list(ct.sum=0, pow.file="")
    
    res <- PS.Main(dat=dat, para=para)
    
    #   res$pval[which(res$pval==0)] = min(res$pval[which(res$pval!=0)],na.rm=T)
    #   res$log.fc[which(res$log.fc == Inf)] = max(res$log.fc[which(!is.infinite(res$log.fc))]) + 1
    #   res$log.fc[which(res$log.fc == -Inf)] = min(res$log.fc[which(!is.infinite(res$log.fc))]) -1
    # res$rnk <- - sign(res$log.fc) * log2(res$pval)
    # res[which(res$rnk == max(res$rnk)),"rnk"] <- res[which(res$rnk == max(res$rnk)),"rnk"] + res[which(res$rnk == max(res$rnk)),"log.fc"]
    res$rnk <- sign(res$log.fc) * res$tt
    rnk <- res[, c("rnk"), drop=FALSE]
    return (list(rnk=rnk, res=res))
}
