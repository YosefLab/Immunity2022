library(openxlsx)
library(edgeR)
library(data.table)
library(parallel)

exp_file <- snakemake@input[["exp"]]
cluster_file <- snakemake@input[["clusters"]]
meta_file <- snakemake@input[["meta"]]
gene_file <- snakemake@input[["gene_symbols"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]

if ("batch_var" %in% names(snakemake@params)){
    batch_var <- snakemake@params[["batch_var"]]
} else {
    batch_var <- NULL
}
print(batch_var)

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

clusters <- read.table(cluster_file, sep="\t", header=TRUE, row.names=1)

meta <- fread(
              paste0('zcat "', meta_file, '"'),
              data.table = FALSE)
rownames(meta) <- meta[, 1]
meta[, 1] <- NULL

counts <- counts[, rownames(clusters)]
meta <- meta[rownames(clusters), , drop = FALSE]

# Filter the genes
valid_genes <- rowSums(counts) > 0
counts <- counts[valid_genes, ]
counts <- data.matrix(counts)


geno_cluster <- as.factor(paste(meta$Genotype, clusters$Cluster, sep = "-"))

de_cluster <- function(cluster_id){

    high_cluster <- paste0("KO-", cluster_id)
    low_cluster <- paste0("ctrl-", cluster_id)

    if (!(high_cluster %in% unique(geno_cluster))){
        return(NULL)
    }

    if (!(low_cluster %in% unique(geno_cluster))){
        return(NULL)
    }

    # Remove cells not in either comparison
    to_keep <- geno_cluster %in% union(low_cluster, high_cluster)
    counts_sub <- counts[, to_keep]
    clusters_sub <- geno_cluster[to_keep]
    meta_sub <- meta[to_keep, ]

    # relabel groups as 'low' and 'high'
    group_labels <- ifelse(clusters_sub %in% low_cluster,
                           "low", "high")
    group_labels <- as.factor(group_labels)
    group_labels <- relevel(group_labels, ref = "low")

    clusters_sub <- data.frame(group = group_labels,
                               row.names = colnames(counts_sub))
    
    if (!is.null(batch_var)){
        batch <- as.factor(meta_sub[[batch_var]])
        clusters_sub$batch <- batch
    }

    y <- DGEList(counts = counts_sub, samples = clusters_sub)
    y <- calcNormFactors(y, method = "none") # just use library size

    # Estimate dispersions
    group <- clusters_sub$group
    if (!is.null(batch_var)) {
        batch <- clusters_sub$batch
        design <- model.matrix(~group + batch)
    } else {
        design <- model.matrix(~group)
    }
    y <- estimateDisp(y, design) # 821 seconds for 7k cells

    fit <- glmFit(y, design)

    lrt <- glmLRT(fit, coef = "grouphigh") # Results stored in lrt$table

    results <- as.data.frame(lrt$table)
    results$FDR <- p.adjust(results$PValue, method = "BH")
    results <- results[order(results$FDR, -1 * abs(results$logFC)), ]

    results$cluster <- cluster_id

    results <- cbind(gene = rownames(results), results)
    rownames(results) <- NULL

    return(results)
}

all_cluster_ids <- unique(clusters$Cluster)

all_cluster_de <- mclapply(all_cluster_ids, de_cluster,
                           mc.cores = min(length(all_cluster_ids), 5)
                          )
all_cluster_de <- all_cluster_de[ !sapply(all_cluster_de, is.null) ]

all_concat <- do.call(rbind, all_cluster_de)

# Augment data with Gene Symbols
if (!is.null(gene_file)){
    symbols <- read.table(gene_file, header = TRUE) # columns geneID, GeneSymbol

    all_concat <- merge(all_concat, symbols, by.x = "gene", by.y = "geneID")

    all_concat <- all_concat[c(8, 1:7)]

    all_concat <- all_concat[order(all_concat$cluster, abs(all_concat$logFC)), ]
}


dir.create(dirname(out_file), showWarnings = FALSE)
write.table(all_concat,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

# Write results to excel also
wb <- createWorkbook()

for (cluster in unique(all_concat$cluster)){
    addWorksheet(wb, cluster)
    sub_data <- all_concat[all_concat$cluster == cluster, ]
    sub_data <- sub_data[order(sub_data$FDR, -1 * abs(sub_data$logFC)), ]
    writeData(wb, sheet = cluster, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
