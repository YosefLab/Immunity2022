library('openxlsx')
source("../shared_scripts/cluster_de.R") # Relative to snakemake file

options(mc.cores = 10)

args <- commandArgs(TRUE)


exp_file <- args[1]
meta_file <- args[2]
qc_file <- args[3]

se <- load_exp_and_clusters(exp_file, meta_file, qc_file)

colData(se)["Cluster"] <- colData(se)["Genotype"]

low_clusters <- "ctrl"
high_clusters <- "KO"

mastResult <- run_mast_de(se, low_clusters, high_clusters,
                          formula = ~group + num_umi)

out_file <- args[4]
out_file_xlsx <- args[5]

out_dir <- dirname(out_file)
dir.create(out_dir, showWarnings = FALSE)

write.table(mastResult,
            file = out_file,
            row.names = FALSE,
            sep = "\t")

write.xlsx(mastResult,
           file = out_file_xlsx,
           sheetName = "KO_vs_ctrl")
