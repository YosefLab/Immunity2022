# This is just an attempt to write an R script file for goseq
# to investigate why the results do not look right
# - One possibliity seems to be that I may have been double-FDR
#   correcting the results in the notebook script (GFP_DE.ipynb

library("goseq")

REPO_BASE <- system2("git", c("rev-parse",  "--show-toplevel"), stdout = TRUE)

tpm_file <- "data_filtered/tpm.txt.gz"
tpm <- read.table(gzfile(tpm_file), sep="\t", header = TRUE, row.names = 1)

# Load gene set
test_file <- "gfp_sets.gmt"
lines <- readLines(test_file)
s_lines <- strsplit(lines, "\t")

gs <- lapply(s_lines, function(x) x[3:length(x)])
names(gs) <- sapply(s_lines, function(x) x[1])

gs <- lapply(gs, function(genes){
    sapply(
        strsplit(genes, ',')
        , function(y) y[1]
    )
})

degenes <- gs[["GFP_th1_ctrl_only"]]

# Gmt-formatted file
sigFile <- paste0(REPO_BASE, "/Signatures/Enrichr/ChEA_2016.txt")
sigFile <- paste0(REPO_BASE, "/Signatures/Enrichr/ARCHS4_TFs_Coexp.txt")
sigFile <- paste0(REPO_BASE, "/Signatures/Enrichr/ENCODE_TF_ChIP-seq_2015.txt")


# Create bias data

bias.data <- log1p(rowMeans(tpm))

gene.vector <- as.integer(rownames(tpm) %in% degenes)
names(gene.vector) <- rownames(tpm)

pwf.counts <- nullp(gene.vector, bias.data = bias.data, plot.fit = TRUE)


# Load the signatures
lines <- readLines(sigFile)
s_lines <- strsplit(lines, "\t")

r <- lapply(s_lines, function(x){
    gsgenes <- sapply(
        strsplit(x[-c(1:2)], ',')
        , function(y) y[1]
    )
    cbind(x[1], gsgenes)
})

msig_sets <- as.data.frame(do.call(rbind, r))
colnames(msig_sets) <- c("Signature", "Genes")

# Run analysis
out <- goseq(pwf.counts, "hg19", "geneSymbol",
    method = "Wallenius",
    #method = "Hypergeometric",
    gene2cat = msig_sets, use_genes_without_cat = TRUE)

rownames(out) <- out$category
out$category <- NULL
out <- out[order(out$over_represented_pvalue), ]
out$FDR <- p.adjust(out$over_represented_pvalue, method = "fdr")

write.table(out, outputFile, sep="\t")
