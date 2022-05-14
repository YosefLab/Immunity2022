# ----------------
# This differs from mast.R in that a multi-level
# model is used
# ----------------

library(jsonlite)
library(SummarizedExperiment)
library(dplyr)
library(data.table)
library(dtplyr) # data.table and dplyr
library(MAST)
library(openxlsx)

#X11.options(type = "Xlib")

if ("filterVar" %in% names(snakemake@params)) {
    filterVar <- snakemake@params[["filterVar"]]
} else {
    filterVar <- NULL
}

if ("filterLevel" %in% names(snakemake@params)) {
    filterLevel <- snakemake@params[["filterLevel"]]
} else {
    filterLevel <- NULL
}

if ("filterVar2" %in% names(snakemake@params)) {
    filterVar2 <- snakemake@params[["filterVar2"]]
} else {
    filterVar2 <- NULL
}

if ("filterLevel2" %in% names(snakemake@params)) {
    filterLevel2 <- snakemake@params[["filterLevel2"]]
} else {
    filterLevel2 <- NULL
}

formula <- as.formula(snakemake@params[["formula"]])

# -----------------------------------------------------------------------------
# Load data!
# -----------------------------------------------------------------------------

tpm_file <- snakemake@input[["tpm"]]
qc_file <- snakemake@input[["qc"]]
meta_file <- snakemake@input[["meta"]]

out_file <- snakemake@output[["out"]]
out_xlsx_file <- snakemake@output[["out_xlsx"]]

tpm <- read.delim(gzfile(tpm_file), row.names = 1)
qc <- read.delim(gzfile(qc_file), row.names = 1)
meta <- read.delim(gzfile(meta_file), row.names = 1)

# Center the qc values
qc <- t(t(qc) - colMeans(qc))

# Re-level meta factors for consistence
meta$GFP <- relevel(meta$GFP, ref = "eGFP_neg")
meta$Genotype <- relevel(meta$Genotype, ref = "ctrl")

# Log the data
logtpm <- data.matrix(log2(tpm + 1))

# Create a summarized experiment
colData <- cbind(qc, meta[rownames(qc), ])


if ("otherMeta" %in% names(snakemake@input)) {
    other_meta_file <- snakemake@input[["otherMeta"]]
    other_meta <- read.delim(gzfile(other_meta_file), row.names = 1)
    colData <- cbind(colData, other_meta[rownames(colData), ])
}


exp <- SummarizedExperiment(assays = list(logtpm = logtpm),
                           colData = colData)

assay(exp) <- assays(exp)$logtpm
metadata(exp)$qc_cols <- colnames(qc)
metadata(exp)$meta_cols <- colnames(meta)


# -----------------------------------------------------------------------------
# Run MAST workflow
# -----------------------------------------------------------------------------

sca <- FromMatrix(assay(exp), colData(exp))

# (Optional) Filter samples

if (!is.null(filterVar)){
    sca <- sca[, (colData(sca)[[filterVar]] == filterLevel)]
}
if (!is.null(filterVar2)){
    sca <- sca[, (colData(sca)[[filterVar2]] == filterLevel2)]
}

options(mc.cores = 10)
zlmCond <- zlm(formula, sca)

# coefficients are listed via colnames(zlmCond@coefC)

options(mc.cores = 1)
stim <- summary(zlmCond, doLRT = TRUE)

summaryTable <- as.data.frame(stim$datatable)

contrasts <- as.character(unique(summaryTable$contrast))
# contrasts <- contrasts[contrasts != "(Intercept)"]
names(contrasts) <- contrasts

groupedResults <- lapply(contrasts,
    function(selectedContrast) {

    tsub <- summaryTable[summaryTable$contrast == selectedContrast, ]
    htsub <- tsub[tsub$component == "H", c("primerid", "Pr(>Chisq)")]
    ctsub <- tsub[tsub$component == "C", c("primerid", "coef", "Pr(>Chisq)")]
    dtsub <- tsub[tsub$component == "D", c("primerid", "coef", "Pr(>Chisq)")]
    lfctsub <- tsub[tsub$component == "logFC", c("primerid", "coef")]

    colnames(htsub) <- c("primerid", "pvalH")
    colnames(ctsub) <- c("primerid", "coefC", "pvalC")
    colnames(dtsub) <- c("primerid", "coefD", "pvalD")
    colnames(lfctsub) <- c("primerid", "logFC")

    htsub["FDR"] <- p.adjust(htsub[["pvalH"]], method = "fdr")

    if (selectedContrast == "(Intercept)"){
        allsub <- merge(ctsub, dtsub)
        allsub$pvalH <- NA
        allsub$FDR <- NA
    } else {
        allsub <- merge(merge(htsub, ctsub), dtsub)
    }

    if (nrow(lfctsub) > 0){
        allsub <- merge(allsub, lfctsub) # No lfc for multiple variable contrasts
    } else {
        allsub$logFC <- 0.0
    }

    allsub["contrast"] <- selectedContrast
    allsub <- allsub[c(
        "primerid", "coefC", "pvalC", "coefD", "pvalD",
        "pvalH", "logFC", "FDR", "contrast"
    )]

    return(allsub)
})

all_concat <- do.call(rbind, groupedResults)

dir.create(dirname(out_file), showWarnings = FALSE)
write.table(all_concat,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

# Write results to excel also
wb <- createWorkbook()

for (contrast in names(groupedResults)){
    sheet <- gsub(":", "X", contrast)
    addWorksheet(wb, sheet)
    sub_data <- groupedResults[[contrast]]
    sub_data <- sub_data[order(sub_data$FDR, -1 * abs(sub_data$coefC)), ]
    writeData(wb, sheet = sheet, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx_file, overwrite = TRUE)
