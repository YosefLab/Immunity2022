library(tidyverse)
library(edgeR)

# Load the data

counts <- read.table("rnaseq/counts.txt", sep="\t", header=TRUE, row.names=1)
meta <- read.table("rnaseq/meta.txt", sep="\t", header=TRUE, row.names=1)

meta['Group'] <- paste(meta$Condition, meta$IL.23R, sep='-')

outliers <- c(
    "S1", "S12", "S19", "S20", "S21", "S33", "S34"
)

y <- DGEList(counts=counts, samples=meta, group=meta$Group)
y <- y[, y$samples$Genotype == 'ctrl']

valid_samples <- setdiff(colnames(y$counts), outliers)

y <- y[, valid_samples]
valid_genes <- filterByExpr(y, min.count = 10)

y <- y[valid_genes, ]

y <- calcNormFactors(y, method='RLE')

# %% Run the Anova-like analysis for the heatmap

# y_sub <- y_sub[, y_sub$samples$Group %in% c("IL-12-negative", "IL-1b+IL-6+IL-23-positive", "IL-1b+IL-6+IL-23-negative")]

y_sub <- y
y_sub$samples$Condition <- dropEmptyLevels(y_sub$samples$Condition)
y_sub$samples$Group <- dropEmptyLevels(y_sub$samples$Group)


design <- model.matrix(~Group, y_sub$samples)


y_sub <- estimateDisp(y_sub, design)

fit <- glmFit(y_sub, design)

lrt <- glmLRT(fit, coef=c(2:7))

results_table <- lrt$table[order(lrt$table$PValue), ]

write.table(results_table, "de_anova.txt", sep="\t", col.names = NA)

norm_factors <- y$samples[, c('lib.size', 'norm.factors')]

write.table(norm_factors, "rle_norm_factors.txt", sep="\t", col.names = NA)


# %% Run a pair-wise analysis (vs 12+21 Th1) to create signatures and evaluate similarity

compare_vs_denom <- function(y, denom){
    y_sub <- y
    y_sub$samples$Condition <- dropEmptyLevels(y_sub$samples$Condition)
    y_sub$samples$Condition <- relevel(y_sub$samples$Condition, denom)

    design <- model.matrix(~Condition, y_sub$samples)
    y_sub <- estimateDisp(y_sub, design)

    fit <- glmFit(y_sub, design)

    coefs <- colnames(design)[2:4]
    results_tables <- lapply(coefs, function(coef) {
        num_comparison <- gsub("Condition", "", coef)
        lrt <- glmLRT(fit, coef = coef)
        results_table <- lrt$table[order(lrt$table$PValue), ]
        results_table["FDR"] <- p.adjust(results_table[["PValue"]], method = "fdr")
        results_table["NumComparison"] <- num_comparison
        results_table["DenomComparison"] <- denom
        results_table["GeneSymbol"] <- rownames(results_table)
        rownames(results_table) <- NULL
        return(results_table)
    })
    return(do.call(rbind, results_tables))
}

results_tables_all <- lapply(
    levels(y$samples$Condition),
    function(cond) {compare_vs_denom(y, cond)}
)

results_tables_all <- do.call(rbind, results_tables_all)

write.table(results_tables_all, "pairwise_de.txt", sep = "\t", row.names = FALSE)
