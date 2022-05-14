source("../../../shared_scripts/cluster_de.R")

exp_file = "../../expression/expression_scaled.txt.gz"
cluster_file = "../clusters.txt"

se = load_exp_and_clusters(exp_file, cluster_file)

middle_group = setdiff(
    se$Cluster %>% unique, 
    c("th17_like", "ccl5_high")
)

# DE for the TH17 group
low_clusters = middle_group
high_clusters = "th17_like"

mastResult = run_mast_de(se, low_clusters, high_clusters)

write.table(mastResult, 
            file=paste0("th17_like_VS_mid_MastResults.txt"),
            row.names=FALSE,
            sep="\t")

# DE for the CCL5 group
low_clusters = middle_group
high_clusters = "ccl5_high"

mastResult = run_mast_de(se, low_clusters, high_clusters)

write.table(mastResult, 
            file=paste0("ccl5_high_VS_mid_MastResults.txt"),
            row.names=FALSE,
            sep="\t")
