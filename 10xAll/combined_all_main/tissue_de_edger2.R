library(openxlsx)
library(edgeR)
library(data.table)
library(parallel)


# exp_file <- "data/expression_counts.feather"
# meta_file <- "data/meta.txt.gz"
# gene_file <- "data/genes.tsv"

exp_file <- snakemake@input[["exp"]]
meta_file <- snakemake@input[["meta"]]
gene_file <- snakemake@input[["gene_symbols"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]
out_cpm <- snakemake@output[["out_cpm"]]

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

meta$Tissue <- relevel(
    as.factor(meta$Tissue),
    ref = "Spleen"
)

meta$Batch <- factor(meta$Batch, levels = c(1, 2))



y <- DGEList(counts = counts)
y <- calcNormFactors(y, method = "none") # just use library size

gene_means <- colMeans(t(counts) / y$samples$lib.size * 1e6)

plot(y$AveLogCPM, log2(gene_means+576.0564))


zz = 2**y$AveLogCPM - gene_means

Batch <- meta$Batch
Tissue <- meta$Tissue
design <- model.matrix(~Batch + Tissue)

options(mc.cores = 1)
# 821 seconds for 7k cells
# 1.37 hours for 30k cells
y <- estimateDisp(y, design)

saveRDS(y, file="edger2_y.rds")

fit <- glmFit(y, design) # 5.2 minutes for 30k cells

coef <- c(
    "TissueColon_IEL",
    "TissueColon_LPL",
    "TissueSI_LPL"
)

# 2.3 minutes for 30k cells
lrt <- glmLRT(fit, coef = coef) # Results in lrt$table

results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method = "BH")
results <- results[order(results$FDR), ]

# results <- cbind(gene = rownames(results), results)
# rownames(results) <- NULL

cpms <- cpmByGroup(y, group=Tissue)

# Augment data with Gene Symbols
if (!is.null(gene_file)){
    symbols <- read.table(gene_file, header = TRUE) # columns geneID, GeneSymbol

    results <- merge(results, symbols, by.x = "row.names", by.y = "geneID")
    rownames(results) <- results$Row.names
    results$Row.names <- NULL

    n <- ncol(results)
    results <- results[c(n, 1:(n - 1))]

    results <- results[order(results$FDR, results$LR * -1), ]

    cpms <- merge(cpms, symbols, by.x = "row.names", by.y = "geneID")
    rownames(cpms) <- cpms$Row.names
    cpms$Row.names <- NULL
    n <- ncol(cpms)
    cpms <- cpms[c(n, 1:(n - 1))]
}

results$Contrast <- "All"

res_concat <- results

dir.create(dirname(out_file), showWarnings = FALSE)

write.table(res_concat,
            file = gzfile(out_file),
            row.names = TRUE,
            col.names = NA,
            sep = "\t")

# Write results to excel also
wb <- createWorkbook()

for (contrast in unique(res_concat$Contrast)){
    addWorksheet(wb, contrast)
    sub_data <- res_concat[res_concat$Contrast == contrast, ]
    sub_data <- sub_data[order(sub_data$FDR, -1 * sub_data$LR), ]
    writeData(wb, sheet = contrast, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)

# Also save the CPM
write.table(cpms,
            file = gzfile(out_cpm),
            row.names = TRUE,
            col.names = NA,
            sep = "\t")
