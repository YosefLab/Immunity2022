# This script runs Tissue DE for the LPL and IEL where only the
# Spleen is used as the denominator

library(openxlsx)
library(edgeR)
library(data.table)
library(parallel)


exp_file <- snakemake@input[["exp"]]
meta_file <- snakemake@input[["meta"]]
gene_file <- snakemake@input[["gene_symbols"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]

# exp_file <- "data/expression_counts.feather"
# meta_file <- "data/meta.txt.gz"
# gene_file <- "data/genes.tsv"
# valid_barcodes_file <- "valid_barcodes.txt"
# 
# out_file <- "TissueDE/tissueDE.txt.gz"
# out_xlsx <- "TissueDE/tissueDE.xlsx"

if (endsWith(exp_file, "txt.gz")){
    counts <- fread(
                    paste0('zcat "', exp_file, '"'),
                    data.table = FALSE)
    rownames(counts) <- counts[, 1]
    counts[, 1] <- NULL
} else if (endsWith(exp_file, "feather")) {
    library(feather)
    counts <- read_feather(exp_file)
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$index
    counts$index <- NULL
}

meta <- fread(
              paste0('zcat "', meta_file, '"'),
              data.table = FALSE)

rownames(meta) <- meta[, 1]
meta[, 1] <- NULL

counts <- counts[, rownames(meta)]

# Filter the genes
valid_genes <- rowSums(counts) > 0
counts <- counts[valid_genes, ]
counts <- data.matrix(counts)


y <- DGEList(counts = counts)
y <- calcNormFactors(y, method = "none") # just use library size


meta$isLPL <- (meta$Tissue == "Colon_LPL") | (meta$Tissue == "SI_LPL")
meta$isIEL <- (meta$Tissue == "Colon_IEL")
meta$isSpleen <- (meta$Tissue == "Spleen")

meta$isLPL <- factor(meta$isLPL, levels = c(FALSE, TRUE))
meta$isIEL <- factor(meta$isIEL, levels = c(FALSE, TRUE))
meta$isSpleen <- factor(meta$isSpleen, levels = c(FALSE, TRUE))

meta$Batch <- factor(meta$Batch, levels = c(1, 2))


options(mc.cores = 1)


## 
Contrasts <- list(LPL = meta$isLPL, IEL = meta$isIEL)
Batch <- meta$Batch
res_all <- lapply(names(Contrasts), function(cc){

    Contrast <- Contrasts[[cc]]

    design <- model.matrix(~Batch + Contrast)

    # Subset design and 'y' to just be contrast + Spleen
    to_keep <- (Contrast == "TRUE") | (meta$isSpleen == "TRUE")
    y <- y[, to_keep]
    design <- design[to_keep, , drop = FALSE]

    # 821 seconds for 7k cells
    # 1.37 hours for 30k cells
    y <- estimateDisp(y, design)

    fit <- glmFit(y, design) # 5.2 minutes for 30k cells

    coef <- "ContrastTRUE"

    # 2.3 minutes for 30k cells
    lrt <- glmLRT(fit, coef = coef) # Results in lrt$table

    results <- as.data.frame(lrt$table)
    results$FDR <- p.adjust(results$PValue, method = "BH")
    results <- results[order(results$FDR), ]

    results <- cbind(gene = rownames(results), results)
    rownames(results) <- NULL

    # Augment data with Gene Symbols
    if (!is.null(gene_file)){
        symbols <- read.table(gene_file, header = TRUE) # columns geneID, GeneSymbol

        results <- merge(results, symbols, by.x = "gene", by.y = "geneID")

        n <- ncol(results)
        results <- results[c(n, 1:(n - 1))]

        results <- results[order(results$FDR, results$LR * -1), ]
    }

    results$Contrast <- cc
    return(results)
})

names(res_all) <- names(Contrasts)

res_concat <- do.call(rbind, res_all)

dir.create(dirname(out_file), showWarnings = FALSE)

write.table(res_concat,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

# Write results to excel also
wb <- createWorkbook()

for (contrast in unique(res_concat$Contrast)){
    addWorksheet(wb, contrast)
    sub_data <- res_concat[res_concat$Contrast == contrast, ]
    sub_data <- sub_data[order(sub_data$FDR, -1 * abs(sub_data$logFC)), ]
    writeData(wb, sheet = contrast, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
