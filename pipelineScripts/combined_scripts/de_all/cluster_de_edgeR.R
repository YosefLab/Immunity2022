library(openxlsx)
library(edgeR)
library(data.table)
library(parallel)

options(mc.cores = 10)

exp_file <- snakemake@input[["exp"]]
meta_file <- snakemake@input[["meta"]]
gene_file <- snakemake@input[["gene_symbols"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]

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

# Estimate dispersions
genotype <- relevel(as.factor(meta$Genotype), ref = "ctrl")
design <- model.matrix(~genotype)
y <- estimateDisp(y, design) # 821 seconds for 7k cells

fit <- glmFit(y, design)

lrt <- glmLRT(fit) # Results in lrt$table

results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method = "BH")
results <- results[order(results$FDR, -1 * abs(results$logFC)), ]

results <- cbind(gene = rownames(results), results)
rownames(results) <- NULL


# Augment data with Gene Symbols
if (!is.null(gene_file)){
    symbols <- read.table(gene_file, header = TRUE) # columns geneID, GeneSymbol

    results <- merge(results, symbols, by.x = "gene", by.y = "geneID")

    results <- results[c(7, 1:6)]

    results <- results[order(abs(results$logFC)), ]
}

results$cluster <- "KO_vs_ctrl"


dir.create(dirname(out_file), showWarnings = FALSE)
write.table(results,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

write.xlsx(results,
           file = out_xlsx,
           sheetName = "KO_vs_ctrl")
