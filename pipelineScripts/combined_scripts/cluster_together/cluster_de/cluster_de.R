library(openxlsx)
library(tidyverse)
library(BiocParallel)
source("../shared_scripts/cluster_de.R")

options(mc.cores = 10)

args <- commandArgs(TRUE)

exp_file <- args[1]
meta_file <- args[2]
qc_file <- args[3]
cluster_file <- args[4]

out_file <- args[5]
out_file_xlsx <- args[6]

se <- load_exp_and_clusters(exp_file, meta_file, qc_file)

clusters <- readr::read_tsv(cluster_file) %>% as.data.frame
rownames(clusters) <- clusters[[1]]
clusters[[1]] <- NULL

colData(se)["CombinedCluster"] <- clusters[rownames(colData(se)), "Cluster"]

colData(se)["Cluster"] <- paste(
                                colData(se)[["Genotype"]],
                                colData(se)[["CombinedCluster"]],
                                sep = "-"
                                )

all_cluster_ids <- unique(colData(se)[["CombinedCluster"]])

de_cluster <- function(cluster){

    high_cluster <- paste0("KO-", cluster)
    low_cluster <- paste0("ctrl-", cluster)

    if (!(high_cluster %in% unique(colData(se))[["Cluster"]])){
        return(NULL)
    }

    if (!(low_cluster %in% unique(colData(se))[["Cluster"]])){
        return(NULL)
    }

    mastResult <- run_mast_de(se, low_cluster, high_cluster,
                              formula = ~group + num_umi)

    mastResult[, "cluster" := rep(cluster, nrow(mastResult))]

    return(mastResult)
}

all_cluster_de <- lapply(all_cluster_ids, de_cluster)

all_cluster_de <- all_cluster_de[ !sapply(all_cluster_de, is.null) ]

all_concat <- do.call(rbind, all_cluster_de)

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

saveWorkbook(wb, file = out_file_xlsx, overwrite = TRUE)
