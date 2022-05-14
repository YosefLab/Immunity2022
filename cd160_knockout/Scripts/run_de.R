library(tidyverse)
library(edgeR)

# Load the data

counts <- read.table("rnaseq/counts.txt", sep="\t", header=TRUE, row.names=1)
meta <- read.table("rnaseq/meta.txt", sep="\t", header=TRUE, row.names=1)

outliers <- "S1"

y <- DGEList(counts=counts, samples=meta, group=meta$Genotype)

valid_samples <- setdiff(colnames(y$counts), outliers)

y <- y[, valid_samples]
valid_genes <- filterByExpr(y, min.count = 10)

y <- y[valid_genes, ]

y <- calcNormFactors(y, method='RLE')

# %% Run the Anova-like analysis for the heatmap

design <- model.matrix(~y$samples$Genotype, y$samples)

y <- estimateDisp(y, design)

fit <- glmFit(y, design)

lrt <- glmLRT(fit, coef=2)

results_table <- lrt$table[order(lrt$table$PValue), ]

results_table['FDR'] <- p.adjust(results_table[['PValue']], method='fdr')

write.table(results_table, "de_genotype.txt", sep="\t", col.names = NA)

norm_factors <- y$samples[, c('lib.size', 'norm.factors')]

write.table(norm_factors, "rle_norm_factors.txt", sep="\t", col.names = NA)
