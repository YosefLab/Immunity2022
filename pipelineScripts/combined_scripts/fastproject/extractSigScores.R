library(feather)
library(VISION)

vis_file <- snakemake@input[["vis"]]
out_file <- snakemake@output[["out"]]

out_dir <- dirname(out_file)
dir.create(out_dir, showWarnings = FALSE)

vis <- readRDS(vis_file)

sigScores <- vis@sigScores
sigScores <- as.data.frame(sigScores)
sigScores$index <- rownames(sigScores)

feather::write_feather(sigScores, out_file)
