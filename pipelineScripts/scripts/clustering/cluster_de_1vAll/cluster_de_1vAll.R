library(openxlsx)
library(dplyr)
library(BiocParallel)
source("../shared_scripts/cluster_de.R")


exp_file <- snakemake@input[["exp"]]
qc_file <- snakemake@input[["qc"]]
cluster_file <- snakemake@input[["clusters"]]

out_file <- snakemake@output[["out"]]
out_xlsx <- snakemake@output[["out_xlsx"]]

se <- load_exp_and_clusters(exp_file, cluster_file, qc_file)


all_cluster_ids <- se$Cluster %>% unique

de_cluster <- function(cluster_id){


    out_group <- setdiff(
        se$Cluster %>% unique,
        cluster_id
    )

    low_clusters <- out_group
    high_clusters <- cluster_id

    mastResult <- run_mast_de_qc(se, low_clusters, high_clusters)

    mastResult[, "cluster" := rep(cluster_id, nrow(mastResult))]

    return(mastResult)
}

param <- BiocParallel::MulticoreParam(
                          workers = min(length(all_cluster_ids), 10)
                          )

all_cluster_de <- bplapply(all_cluster_ids, de_cluster,
                           BPPARAM = param)

all_concat <- do.call(rbind, all_cluster_de)

dir.create(dirname(out_file), showWarnings = FALSE)
write.table(all_concat,
            file = gzfile(out_file),
            row.names = FALSE,
            sep = "\t")

# Write results to excel also
all_concat <- as.data.frame(all_concat)
wb <- createWorkbook()

for (cluster in unique(all_concat$cluster)){
    addWorksheet(wb, cluster)
    sub_data <- all_concat[all_concat$cluster == cluster, ]
    writeData(wb, sheet = cluster, x = sub_data)
}

saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
