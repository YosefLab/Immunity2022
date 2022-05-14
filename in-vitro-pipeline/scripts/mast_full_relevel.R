# -------------
# This is the same as mast_full except 'ctrl = 1' in the Genotype encoding
#
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

# Re-level meta factors for consistency
meta$GFP <- relevel(meta$GFP, ref = "eGFP_neg")
meta$Genotype <- relevel(meta$Genotype, ref = "KO")

# Log the data
logtpm <- data.matrix(log2(tpm + 1))

# Create a summarized experiment
colData <- cbind(qc, meta[rownames(qc), ])


if ("otherMeta" %in% names(snakemake@input)) {
    other_meta_file <- snakemake@input[["otherMeta"]]
    other_meta <- read.delim(gzfile(other_meta_file), row.names = 1)
    # Center other-meta too
    other_meta <- t(t(other_meta) - colMeans(other_meta))
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

doLRT <- list(
    "Genotypectrl:GFPeGFP_pos",
    c("GFPeGFP_pos", "Genotypectrl:GFPeGFP_pos"),
    "GFPeGFP_pos"
    )

stim <- summary(zlmCond, doLRT = doLRT)

summaryTable <- as.data.frame(stim$datatable)


levels(summaryTable$contrast)[
    levels(summaryTable$contrast) == "c(\"GFPeGFP_pos\", \"Genotypectrl:GFPeGFP_pos\")"
    ] <- "GFP+Genotype:GFP"


contrasts <- c("GFPeGFP_pos", "GFP+Genotype:GFP", "Genotypectrl:GFPeGFP_pos")

groupedResults <- lapply(setNames(contrasts, contrasts),
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

    allsub <- merge(merge(htsub, ctsub), dtsub)
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

# Now, the logFC are condition-dependent.  So we need to fix them by
# evaluating it in specific conditions.  Do this with two conditions:
#    GFP_pos, GenotypeKO
#    GFP_pos, GenotypeCtrl
# Omitting GFP_neg as it's not relevant

coefnames <- colnames(coef(zlmCond, 'D'))
contrast0 <- setNames(rep(0, length(coefnames)), coefnames)
contrast0[c("(Intercept)", "GFPeGFP_pos")] <- 1

contrast1 <- contrast0
contrast1[c("Genotypectrl:GFPeGFP_pos")] <- 1

lfc_ix_ko <- logFC(zlmCond, contrast0, contrast1)
lfc_ix_ko <- lfc_ix_ko$logFC
lfc_ix_ko <- setNames(lfc_ix_ko[, 1], rownames(lfc_ix_ko))

contrast0[c("Genotypectrl")] <- 1
contrast1[c("Genotypectrl")] <- 1

lfc_ix_ctrl <- logFC(zlmCond, contrast0, contrast1)
lfc_ix_ctrl <- lfc_ix_ctrl$logFC
lfc_ix_ctrl <- setNames(lfc_ix_ctrl[, 1], rownames(lfc_ix_ctrl))

lfc_ix <- lfc_ix_ko / 2 + lfc_ix_ctrl / 2

ix_results <- groupedResults[["Genotypectrl:GFPeGFP_pos"]]
lfc_ix <- lfc_ix[ix_results$primerid]
ix_results$logFC <- lfc_ix
groupedResults[["Genotypectrl:GFPeGFP_pos"]] <- ix_results


# Now do the heirchical selection
GFP_res <- groupedResults[["GFP+Genotype:GFP"]]
GFP_only_res <- groupedResults[["GFPeGFP_pos"]]
sig_gfp <- GFP_res[GFP_res$FDR < 0.1, "primerid"]

interaction <- groupedResults[["Genotypectrl:GFPeGFP_pos"]]
interaction <- interaction[
    interaction$primerid %in% sig_gfp,
    ]
interaction["FDR"] <- p.adjust(interaction[["pvalH"]], method = "fdr")

sig_dep <- interaction[interaction$FDR < 0.1, "primerid"]

sig_ind <- setdiff(sig_gfp, sig_dep)
GFP_only_ind <- GFP_only_res[
    GFP_only_res$primerid %in% sig_ind, ]
GFP_only_ind$contrast <- "GFP_Ind"

groupedResults <- list(
    "GFPeGFP_pos" = GFP_only_res,
    "GFP+Genotype:GFP" = GFP_res,
    "Genotypectrl:GFPeGFP_pos" = groupedResults[["Genotypectrl:GFPeGFP_pos"]],
    "Genotypectrl:GFPeGFP_pos_sub" = interaction,
    "GFP_Ind_sub" = GFP_only_ind
    )

groupedResults_x <- lapply(names(groupedResults), function(contrast){
    res <- groupedResults[[contrast]]
    res$contrast <- contrast
    return(res)
})

names(groupedResults_x) <- names(groupedResults)
groupedResults <- groupedResults_x

# Results in the following tabs
#            GFPeGFP_pos:  1. Test GFP coef
#       GFP+Genotype:GFP:  2. Test GFP & Interaction coef
# Genotypectrl:GFPeGFP_pos:  3. Test interaction coef, subset for sig in #2
#                GFP_Ind:  4. Test GFP coef, subset for #2 - #3


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
    sub_data <- sub_data[order(sub_data$FDR), ]
    writeData(wb, sheet = sheet, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx_file, overwrite = TRUE)
