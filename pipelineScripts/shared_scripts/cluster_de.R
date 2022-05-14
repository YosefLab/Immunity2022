library(tidyverse) # ggplot2, tibble, tidyr, readr, purrr, dplyr
library(SummarizedExperiment)
library(MAST)
library(data.table)
library(dtplyr)
library(feather)
library(tools)

load_exp_and_clusters <- function(exp_file, meta_file, qc_file=NULL){

    if (tools::file_ext(exp_file) == "feather") {
        exp <- feather::read_feather(exp_file)
        exp <- as.data.frame(exp)
        rownames(exp) <- exp$index
        exp$index <- NULL

    } else {
        exp <- readr::read_tsv(exp_file) %>% as.data.frame
        rownames(exp) <- exp$GeneSymbol
        exp$GeneSymbol <- NULL
    }


    exp <- data.matrix(exp)

    clusters = readr::read_tsv(meta_file) %>% as.data.frame
    rownames(clusters) = clusters[[1]]
    clusters[[1]] = NULL

    # Log-transform data
    logexp = log2(exp+1)

    # align clusters with log data
    clusters = clusters[colnames(logexp), , drop=FALSE]

    if(!is.null(qc_file)){
        qc = readr::read_tsv(qc_file) %>% as.data.frame
        rownames(qc) = qc[[1]]
        qc[[1]]= NULL

        colData = merge(qc, clusters, by="row.names", all=TRUE)
        rownames(colData) = colData$Row.names
        colData$Row.names = NULL

        # align colData with log data
        colData = colData[colnames(logexp), , drop=FALSE]
    } else {
        colData = clusters
    }

    se = SummarizedExperiment(assays=list(logcpm=as.matrix(logexp)),
                              colData=colData)

    return(se)

}


run_mast_de = function(se, low_clusters, high_clusters, formula=~group){

    # subset data for only those in each group
    to_keep = se$Cluster %in% union(low_clusters, high_clusters)

    se_sub = se[, to_keep]

    # relabel groups as 'low' and 'high'
    group_labels = ifelse(se_sub$Cluster %in% low_clusters, "low", "high")
    group_labels = group_labels %>% as.factor %>% relevel(ref="low")

    se_sub@colData$group = group_labels

    # define MAST SCA object
    sca = MAST::FromMatrix(assay(se_sub), colData(se_sub))

    # run test
    zlmCond = zlm(formula, sca)
    stim = summary(zlmCond, doLRT="grouphigh")

    # transform results
    toSigTable = function(summaryCond, selectedContrast){
        summaryDt <- summaryCond$datatable
        fcHurdle <- merge(
                      summaryDt[contrast==selectedContrast & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==selectedContrast & component=='logFC',
                                .(primerid, coef, ci.hi, ci.lo)]
                      , by='primerid') #logFC coefficients

        dcoef = summaryDt[contrast==selectedContrast & component=='D',
                  .(primerid, coef)] #hurdle P values

        colnames(dcoef) = c("primerid", "dcoef")

        fcHurdle = merge(fcHurdle, dcoef, by='primerid')

        fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
        setorder(fcHurdle, fdr)

        # set a proxy-coef which is coefC unless coefC is nan, in which
        # case it's coefD

        fcHurdle$proxyCoef = fcHurdle[,'coef']
        na_proxy = fcHurdle[is.na(proxyCoef)]$primerid
        fcHurdle[is.na(proxyCoef), proxyCoef := dcoef]

        setnames(fcHurdle, old='Pr(>Chisq)', new='pval')

        return(fcHurdle)
    }

    stimSig = toSigTable(stim, "grouphigh")

    return(stimSig)
}


# Same as above, but use qPC1 and qPC2 as part of the formula
run_mast_de_qc = function(se, low_clusters, high_clusters){

    # subset data for only those in each group
    to_keep = se$Cluster %in% union(low_clusters, high_clusters)

    se_sub = se[, to_keep]

    # relabel groups as 'low' and 'high'
    group_labels = ifelse(se_sub$Cluster %in% low_clusters, "low", "high")
    group_labels = group_labels %>% as.factor %>% relevel(ref="low")

    se_sub@colData$group = group_labels

    # define MAST SCA object
    sca = MAST::FromMatrix(assay(se_sub), colData(se_sub))

    # run test
    zlmCond = zlm(~group + qPC1 + qPC2, sca)
    stim = summary(zlmCond, doLRT="grouphigh")

    # transform results
    toSigTable = function(summaryCond, selectedContrast){
        summaryDt <- summaryCond$datatable
        fcHurdle <- merge(
                      summaryDt[contrast==selectedContrast & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==selectedContrast & component=='logFC',
                                .(primerid, coef, ci.hi, ci.lo)]
                      , by='primerid') #logFC coefficients

        dcoef = summaryDt[contrast==selectedContrast & component=='D',
                  .(primerid, coef)] #hurdle P values

        colnames(dcoef) = c("primerid", "dcoef")

        fcHurdle = merge(fcHurdle, dcoef, by='primerid')

        fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
        setorder(fcHurdle, fdr)

        # set a proxy-coef which is coefC unless coefC is nan, in which
        # case it's coefD

        fcHurdle$proxyCoef = fcHurdle[,'coef']
        na_proxy = fcHurdle[is.na(proxyCoef)]$primerid
        fcHurdle[is.na(proxyCoef), proxyCoef := dcoef]

        setnames(fcHurdle, old='Pr(>Chisq)', new='pval')

        return(fcHurdle)
    }

    stimSig = toSigTable(stim, "grouphigh")

    return(stimSig)
}
