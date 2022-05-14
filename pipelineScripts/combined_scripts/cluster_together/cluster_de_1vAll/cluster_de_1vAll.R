library(openxlsx)
library(dplyr)
library(BiocParallel)
source("../../../shared_scripts/cluster_de.R")

options(mc.cores = 10)

exp_file <- "../../data/expression_scaled.txt.gz"
cluster_file <- "../clusters.txt"
qc_file <- "../../data/qc.txt.gz"

se <- load_exp_and_clusters(exp_file, cluster_file, qc_file)

all_cluster_ids <- se$Cluster %>% unique

de_cluster <- function(cluster_id){


    out_group <- setdiff(
        se$Cluster %>% unique, 
        cluster_id
    )

    low_clusters <- out_group
    high_clusters <- cluster_id

    mastResult <- run_mast_de(se, low_clusters, high_clusters, formula=~group + num_umi)

    mastResult[, "cluster" := rep(cluster_id, nrow(mastResult))]

    return(mastResult)
}

all_cluster_de <- lapply(all_cluster_ids, de_cluster)


all_concat <- do.call(rbind, all_cluster_de)

write.table(all_concat,
            file = gzfile("cluster_de_1vAll_Results.txt.gz"),
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

saveWorkbook(wb, file = "cluster_de_1vAll_Results.xlsx", overwrite = TRUE)
